// Provides an efficient blocking version of mcl::concurrent_queue.
// ©2015-2020 Cameron Desrochers. Distributed under the terms of the simplified
// BSD license, available at the top of concurrent_queue.h.
// Also dual-licensed under the Boost Software License (see LICENSE.md)
// Uses Jeff Preshing's semaphore implementation (under the terms of its
// separate zlib license, see lightweight_semaphore.h).

#pragma once

#include "concurrent_queue.h"
#include "lightweight_semaphore.h"

#include <type_traits>
#include <cerrno>
#include <memory>
#include <chrono>
#include <ctime>

namespace mcl
{
// This is a blocking version of the queue. It has an almost identical interface to
// the normal non-blocking version, with the addition of various wait_pop() methods
// and the removal of producer-specific dequeue methods.
template<typename T, typename Traits = ConcurrentQueueDefaultTraits>
class blocking_concurrent_queue
{
private:
	typedef ::mcl::concurrent_queue<T, Traits> concurrent_queue;
	typedef ::mcl::lightweight_semaphore lightweight_semaphore;

public:
	typedef typename concurrent_queue::producer_token_t producer_token_t;
	typedef typename concurrent_queue::consumer_token_t consumer_token_t;
	
	typedef typename concurrent_queue::index_t index_t;
	typedef typename concurrent_queue::size_t size_t;
	typedef typename std::make_signed<size_t>::type ssize_t;
	
	static const size_t BLOCK_SIZE = concurrent_queue::BLOCK_SIZE;
	static const size_t EXPLICIT_BLOCK_EMPTY_COUNTER_THRESHOLD = concurrent_queue::EXPLICIT_BLOCK_EMPTY_COUNTER_THRESHOLD;
	static const size_t EXPLICIT_INITIAL_INDEX_SIZE = concurrent_queue::EXPLICIT_INITIAL_INDEX_SIZE;
	static const size_t IMPLICIT_INITIAL_INDEX_SIZE = concurrent_queue::IMPLICIT_INITIAL_INDEX_SIZE;
	static const size_t INITIAL_IMPLICIT_PRODUCER_HASH_SIZE = concurrent_queue::INITIAL_IMPLICIT_PRODUCER_HASH_SIZE;
	static const std::uint32_t EXPLICIT_CONSUMER_CONSUMPTION_QUOTA_BEFORE_ROTATE = concurrent_queue::EXPLICIT_CONSUMER_CONSUMPTION_QUOTA_BEFORE_ROTATE;
	static const size_t MAX_SUBQUEUE_SIZE = concurrent_queue::MAX_SUBQUEUE_SIZE;
	
public:
	// Creates a queue with at least `capacity` element slots; note that the
	// actual number of elements that can be inserted without additional memory
	// allocation depends on the number of producers and the block size (e.g. if
	// the block size is equal to `capacity`, only a single block will be allocated
	// up-front, which means only a single producer will be able to enqueue elements
	// without an extra allocation -- blocks aren't shared between producers).
	// This method is not thread safe -- it is up to the user to ensure that the
	// queue is fully constructed before it starts being used by other threads (this
	// includes making the memory effects of construction visible, possibly with a
	// memory barrier).
	explicit blocking_concurrent_queue(size_t capacity = 6 * BLOCK_SIZE)
		: inner(capacity), sema(create<lightweight_semaphore, ssize_t, int>(0, (int)Traits::MAX_SEMA_SPINS), &blocking_concurrent_queue::template destroy<lightweight_semaphore>)
	{
		assert(reinterpret_cast<concurrent_queue*>((blocking_concurrent_queue*)1) == &((blocking_concurrent_queue*)1)->inner && "blocking_concurrent_queue must have concurrent_queue as its first member");
		if (!sema) {
			MOODYCAMEL_THROW(std::bad_alloc());
		}
	}
	
	blocking_concurrent_queue(size_t minCapacity, size_t maxExplicitProducers, size_t maxImplicitProducers)
		: inner(minCapacity, maxExplicitProducers, maxImplicitProducers), sema(create<lightweight_semaphore, ssize_t, int>(0, (int)Traits::MAX_SEMA_SPINS), &blocking_concurrent_queue::template destroy<lightweight_semaphore>)
	{
		assert(reinterpret_cast<concurrent_queue*>((blocking_concurrent_queue*)1) == &((blocking_concurrent_queue*)1)->inner && "blocking_concurrent_queue must have concurrent_queue as its first member");
		if (!sema) {
			MOODYCAMEL_THROW(std::bad_alloc());
		}
	}
	
	// Disable copying and copy assignment
	blocking_concurrent_queue(blocking_concurrent_queue const&) MOODYCAMEL_DELETE_FUNCTION;
	blocking_concurrent_queue& operator=(blocking_concurrent_queue const&) MOODYCAMEL_DELETE_FUNCTION;
	
	// Moving is supported, but note that it is *not* a thread-safe operation.
	// Nobody can use the queue while it's being moved, and the memory effects
	// of that move must be propagated to other threads before they can use it.
	// Note: When a queue is moved, its tokens are still valid but can only be
	// used with the destination queue (i.e. semantically they are moved along
	// with the queue itself).
	blocking_concurrent_queue(blocking_concurrent_queue&& other) MOODYCAMEL_NOEXCEPT
		: inner(std::move(other.inner)), sema(std::move(other.sema))
	{ }
	
	inline blocking_concurrent_queue& operator=(blocking_concurrent_queue&& other) MOODYCAMEL_NOEXCEPT
	{
		return swap_internal(other);
	}
	
	// Swaps this queue's state with the other's. Not thread-safe.
	// Swapping two queues does not invalidate their tokens, however
	// the tokens that were created for one queue must be used with
	// only the swapped queue (i.e. the tokens are tied to the
	// queue's movable state, not the object itself).
	inline void swap(blocking_concurrent_queue& other) MOODYCAMEL_NOEXCEPT
	{
		swap_internal(other);
	}
	
private:
	blocking_concurrent_queue& swap_internal(blocking_concurrent_queue& other)
	{
		if (this == &other) {
			return *this;
		}
		
		inner.swap(other.inner);
		sema.swap(other.sema);
		return *this;
	}
	
public:
	// Enqueues a single item (by copying it).
	// Allocates memory if required. Only fails if memory allocation fails (or implicit
	// production is disabled because Traits::INITIAL_IMPLICIT_PRODUCER_HASH_SIZE is 0,
	// or Traits::MAX_SUBQUEUE_SIZE has been defined and would be surpassed).
	// Thread-safe.
	inline bool push(T const& item)
	{
		if ((details::likely)(inner.push(item))) {
			sema->signal();
			return true;
		}
		return false;
	}
	
	// Enqueues a single item (by moving it, if possible).
	// Allocates memory if required. Only fails if memory allocation fails (or implicit
	// production is disabled because Traits::INITIAL_IMPLICIT_PRODUCER_HASH_SIZE is 0,
	// or Traits::MAX_SUBQUEUE_SIZE has been defined and would be surpassed).
	// Thread-safe.
	inline bool push(T&& item)
	{
		if ((details::likely)(inner.push(std::move(item)))) {
			sema->signal();
			return true;
		}
		return false;
	}
	
	// Enqueues a single item (by copying it) using an explicit producer token.
	// Allocates memory if required. Only fails if memory allocation fails (or
	// Traits::MAX_SUBQUEUE_SIZE has been defined and would be surpassed).
	// Thread-safe.
	inline bool push(producer_token_t const& token, T const& item)
	{
		if ((details::likely)(inner.push(token, item))) {
			sema->signal();
			return true;
		}
		return false;
	}
	
	// Enqueues a single item (by moving it, if possible) using an explicit producer token.
	// Allocates memory if required. Only fails if memory allocation fails (or
	// Traits::MAX_SUBQUEUE_SIZE has been defined and would be surpassed).
	// Thread-safe.
	inline bool push(producer_token_t const& token, T&& item)
	{
		if ((details::likely)(inner.push(token, std::move(item)))) {
			sema->signal();
			return true;
		}
		return false;
	}
	
	// Enqueues several items.
	// Allocates memory if required. Only fails if memory allocation fails (or
	// implicit production is disabled because Traits::INITIAL_IMPLICIT_PRODUCER_HASH_SIZE
	// is 0, or Traits::MAX_SUBQUEUE_SIZE has been defined and would be surpassed).
	// Note: Use std::make_move_iterator if the elements should be moved instead of copied.
	// Thread-safe.
	template<typename It>
	inline bool push_bulk(It itemFirst, size_t count)
	{
		if ((details::likely)(inner.push_bulk(std::forward<It>(itemFirst), count))) {
			sema->signal((lightweight_semaphore::ssize_t)(ssize_t)count);
			return true;
		}
		return false;
	}
	
	// Enqueues several items using an explicit producer token.
	// Allocates memory if required. Only fails if memory allocation fails
	// (or Traits::MAX_SUBQUEUE_SIZE has been defined and would be surpassed).
	// Note: Use std::make_move_iterator if the elements should be moved
	// instead of copied.
	// Thread-safe.
	template<typename It>
	inline bool push_bulk(producer_token_t const& token, It itemFirst, size_t count)
	{
		if ((details::likely)(inner.push_bulk(token, std::forward<It>(itemFirst), count))) {
			sema->signal((lightweight_semaphore::ssize_t)(ssize_t)count);
			return true;
		}
		return false;
	}
	
	// Enqueues a single item (by copying it).
	// Does not allocate memory. Fails if not enough room to enqueue (or implicit
	// production is disabled because Traits::INITIAL_IMPLICIT_PRODUCER_HASH_SIZE
	// is 0).
	// Thread-safe.
	inline bool try_push(T const& item)
	{
		if (inner.try_push(item)) {
			sema->signal();
			return true;
		}
		return false;
	}
	
	// Enqueues a single item (by moving it, if possible).
	// Does not allocate memory (except for one-time implicit producer).
	// Fails if not enough room to enqueue (or implicit production is
	// disabled because Traits::INITIAL_IMPLICIT_PRODUCER_HASH_SIZE is 0).
	// Thread-safe.
	inline bool try_push(T&& item)
	{
		if (inner.try_push(std::move(item))) {
			sema->signal();
			return true;
		}
		return false;
	}
	
	// Enqueues a single item (by copying it) using an explicit producer token.
	// Does not allocate memory. Fails if not enough room to enqueue.
	// Thread-safe.
	inline bool try_push(producer_token_t const& token, T const& item)
	{
		if (inner.try_push(token, item)) {
			sema->signal();
			return true;
		}
		return false;
	}
	
	// Enqueues a single item (by moving it, if possible) using an explicit producer token.
	// Does not allocate memory. Fails if not enough room to enqueue.
	// Thread-safe.
	inline bool try_push(producer_token_t const& token, T&& item)
	{
		if (inner.try_push(token, std::move(item))) {
			sema->signal();
			return true;
		}
		return false;
	}
	
	// Enqueues several items.
	// Does not allocate memory (except for one-time implicit producer).
	// Fails if not enough room to enqueue (or implicit production is
	// disabled because Traits::INITIAL_IMPLICIT_PRODUCER_HASH_SIZE is 0).
	// Note: Use std::make_move_iterator if the elements should be moved
	// instead of copied.
	// Thread-safe.
	template<typename It>
	inline bool try_push_bulk(It itemFirst, size_t count)
	{
		if (inner.try_push_bulk(std::forward<It>(itemFirst), count)) {
			sema->signal((lightweight_semaphore::ssize_t)(ssize_t)count);
			return true;
		}
		return false;
	}
	
	// Enqueues several items using an explicit producer token.
	// Does not allocate memory. Fails if not enough room to enqueue.
	// Note: Use std::make_move_iterator if the elements should be moved
	// instead of copied.
	// Thread-safe.
	template<typename It>
	inline bool try_push_bulk(producer_token_t const& token, It itemFirst, size_t count)
	{
		if (inner.try_push_bulk(token, std::forward<It>(itemFirst), count)) {
			sema->signal((lightweight_semaphore::ssize_t)(ssize_t)count);
			return true;
		}
		return false;
	}
	
	
	// Attempts to dequeue from the queue.
	// Returns false if all producer streams appeared empty at the time they
	// were checked (so, the queue is likely but not guaranteed to be empty).
	// Never allocates. Thread-safe.
	template<typename U>
	inline bool try_pop(U& item)
	{
		if (sema->tryWait()) {
			while (!inner.try_pop(item)) {
				continue;
			}
			return true;
		}
		return false;
	}
	
	// Attempts to dequeue from the queue using an explicit consumer token.
	// Returns false if all producer streams appeared empty at the time they
	// were checked (so, the queue is likely but not guaranteed to be empty).
	// Never allocates. Thread-safe.
	template<typename U>
	inline bool try_pop(consumer_token_t& token, U& item)
	{
		if (sema->tryWait()) {
			while (!inner.try_pop(token, item)) {
				continue;
			}
			return true;
		}
		return false;
	}
	
	// Attempts to dequeue several elements from the queue.
	// Returns the number of items actually dequeued.
	// Returns 0 if all producer streams appeared empty at the time they
	// were checked (so, the queue is likely but not guaranteed to be empty).
	// Never allocates. Thread-safe.
	template<typename It>
	inline size_t try_pop_bulk(It itemFirst, size_t max)
	{
		size_t count = 0;
		max = (size_t)sema->tryWaitMany((lightweight_semaphore::ssize_t)(ssize_t)max);
		while (count != max) {
			count += inner.template try_pop_bulk<It&>(itemFirst, max - count);
		}
		return count;
	}
	
	// Attempts to dequeue several elements from the queue using an explicit consumer token.
	// Returns the number of items actually dequeued.
	// Returns 0 if all producer streams appeared empty at the time they
	// were checked (so, the queue is likely but not guaranteed to be empty).
	// Never allocates. Thread-safe.
	template<typename It>
	inline size_t try_pop_bulk(consumer_token_t& token, It itemFirst, size_t max)
	{
		size_t count = 0;
		max = (size_t)sema->tryWaitMany((lightweight_semaphore::ssize_t)(ssize_t)max);
		while (count != max) {
			count += inner.template try_pop_bulk<It&>(token, itemFirst, max - count);
		}
		return count;
	}
	
	
	
	// Blocks the current thread until there's something to dequeue, then
	// dequeues it.
	// Never allocates. Thread-safe.
	template<typename U>
	inline void wait_pop(U& item)
	{
		while (!sema->wait()) {
			continue;
		}
		while (!inner.try_pop(item)) {
			continue;
		}
	}

	// Blocks the current thread until either there's something to dequeue
	// or the timeout (specified in microseconds) expires. Returns false
	// without setting `item` if the timeout expires, otherwise assigns
	// to `item` and returns true.
	// Using a negative timeout indicates an indefinite timeout,
	// and is thus functionally equivalent to calling wait_pop.
	// Never allocates. Thread-safe.
	template<typename U>
	inline bool wait_pop_timed(U& item, std::int64_t timeout_usecs)
	{
		if (!sema->wait(timeout_usecs)) {
			return false;
		}
		while (!inner.try_pop(item)) {
			continue;
		}
		return true;
	}
    
    // Blocks the current thread until either there's something to dequeue
	// or the timeout expires. Returns false without setting `item` if the
    // timeout expires, otherwise assigns to `item` and returns true.
	// Never allocates. Thread-safe.
	template<typename U, typename Rep, typename Period>
	inline bool wait_pop_timed(U& item, std::chrono::duration<Rep, Period> const& timeout)
    {
        return wait_pop_timed(item, std::chrono::duration_cast<std::chrono::microseconds>(timeout).count());
    }
	
	// Blocks the current thread until there's something to dequeue, then
	// dequeues it using an explicit consumer token.
	// Never allocates. Thread-safe.
	template<typename U>
	inline void wait_pop(consumer_token_t& token, U& item)
	{
		while (!sema->wait()) {
			continue;
		}
		while (!inner.try_pop(token, item)) {
			continue;
		}
	}
	
	// Blocks the current thread until either there's something to dequeue
	// or the timeout (specified in microseconds) expires. Returns false
	// without setting `item` if the timeout expires, otherwise assigns
	// to `item` and returns true.
	// Using a negative timeout indicates an indefinite timeout,
	// and is thus functionally equivalent to calling wait_pop.
	// Never allocates. Thread-safe.
	template<typename U>
	inline bool wait_pop_timed(consumer_token_t& token, U& item, std::int64_t timeout_usecs)
	{
		if (!sema->wait(timeout_usecs)) {
			return false;
		}
		while (!inner.try_pop(token, item)) {
			continue;
		}
		return true;
	}
    
    // Blocks the current thread until either there's something to dequeue
	// or the timeout expires. Returns false without setting `item` if the
    // timeout expires, otherwise assigns to `item` and returns true.
	// Never allocates. Thread-safe.
	template<typename U, typename Rep, typename Period>
	inline bool wait_pop_timed(consumer_token_t& token, U& item, std::chrono::duration<Rep, Period> const& timeout)
    {
        return wait_pop_timed(token, item, std::chrono::duration_cast<std::chrono::microseconds>(timeout).count());
    }
	
	// Attempts to dequeue several elements from the queue.
	// Returns the number of items actually dequeued, which will
	// always be at least one (this method blocks until the queue
	// is non-empty) and at most max.
	// Never allocates. Thread-safe.
	template<typename It>
	inline size_t wait_pop_bulk(It itemFirst, size_t max)
	{
		size_t count = 0;
		max = (size_t)sema->waitMany((lightweight_semaphore::ssize_t)(ssize_t)max);
		while (count != max) {
			count += inner.template try_pop_bulk<It&>(itemFirst, max - count);
		}
		return count;
	}
	
	// Attempts to dequeue several elements from the queue.
	// Returns the number of items actually dequeued, which can
	// be 0 if the timeout expires while waiting for elements,
	// and at most max.
	// Using a negative timeout indicates an indefinite timeout,
	// and is thus functionally equivalent to calling wait_pop_bulk.
	// Never allocates. Thread-safe.
	template<typename It>
	inline size_t wait_pop_bulk_timed(It itemFirst, size_t max, std::int64_t timeout_usecs)
	{
		size_t count = 0;
		max = (size_t)sema->waitMany((lightweight_semaphore::ssize_t)(ssize_t)max, timeout_usecs);
		while (count != max) {
			count += inner.template try_pop_bulk<It&>(itemFirst, max - count);
		}
		return count;
	}
    
    // Attempts to dequeue several elements from the queue.
	// Returns the number of items actually dequeued, which can
	// be 0 if the timeout expires while waiting for elements,
	// and at most max.
	// Never allocates. Thread-safe.
	template<typename It, typename Rep, typename Period>
	inline size_t wait_pop_bulk_timed(It itemFirst, size_t max, std::chrono::duration<Rep, Period> const& timeout)
    {
        return wait_pop_bulk_timed<It&>(itemFirst, max, std::chrono::duration_cast<std::chrono::microseconds>(timeout).count());
    }
	
	// Attempts to dequeue several elements from the queue using an explicit consumer token.
	// Returns the number of items actually dequeued, which will
	// always be at least one (this method blocks until the queue
	// is non-empty) and at most max.
	// Never allocates. Thread-safe.
	template<typename It>
	inline size_t wait_pop_bulk(consumer_token_t& token, It itemFirst, size_t max)
	{
		size_t count = 0;
		max = (size_t)sema->waitMany((lightweight_semaphore::ssize_t)(ssize_t)max);
		while (count != max) {
			count += inner.template try_pop_bulk<It&>(token, itemFirst, max - count);
		}
		return count;
	}
	
	// Attempts to dequeue several elements from the queue using an explicit consumer token.
	// Returns the number of items actually dequeued, which can
	// be 0 if the timeout expires while waiting for elements,
	// and at most max.
	// Using a negative timeout indicates an indefinite timeout,
	// and is thus functionally equivalent to calling wait_pop_bulk.
	// Never allocates. Thread-safe.
	template<typename It>
	inline size_t wait_pop_bulk_timed(consumer_token_t& token, It itemFirst, size_t max, std::int64_t timeout_usecs)
	{
		size_t count = 0;
		max = (size_t)sema->waitMany((lightweight_semaphore::ssize_t)(ssize_t)max, timeout_usecs);
		while (count != max) {
			count += inner.template try_pop_bulk<It&>(token, itemFirst, max - count);
		}
		return count;
	}
	
	// Attempts to dequeue several elements from the queue using an explicit consumer token.
	// Returns the number of items actually dequeued, which can
	// be 0 if the timeout expires while waiting for elements,
	// and at most max.
	// Never allocates. Thread-safe.
	template<typename It, typename Rep, typename Period>
	inline size_t wait_pop_bulk_timed(consumer_token_t& token, It itemFirst, size_t max, std::chrono::duration<Rep, Period> const& timeout)
    {
        return wait_pop_bulk_timed<It&>(token, itemFirst, max, std::chrono::duration_cast<std::chrono::microseconds>(timeout).count());
    }
	
	
	// Returns an estimate of the total number of elements currently in the queue. This
	// estimate is only accurate if the queue has completely stabilized before it is called
	// (i.e. all enqueue and dequeue operations have completed and their memory effects are
	// visible on the calling thread, and no further operations start while this method is
	// being called).
	// Thread-safe.
	inline size_t size_approx() const
	{
		return (size_t)sema->availableApprox();
	}
	
	
	// Returns true if the underlying atomic variables used by
	// the queue are lock-free (they should be on most platforms).
	// Thread-safe.
	static constexpr bool is_lock_free()
	{
		return concurrent_queue::is_lock_free();
	}
	

private:
	template<typename U, typename A1, typename A2>
	static inline U* create(A1&& a1, A2&& a2)
	{
		void* p = (Traits::malloc)(sizeof(U));
		return p != nullptr ? new (p) U(std::forward<A1>(a1), std::forward<A2>(a2)) : nullptr;
	}
	
	template<typename U>
	static inline void destroy(U* p)
	{
		if (p != nullptr) {
			p->~U();
		}
		(Traits::free)(p);
	}
	
private:
	concurrent_queue inner;
	std::unique_ptr<lightweight_semaphore, void (*)(lightweight_semaphore*)> sema;
};


template<typename T, typename Traits>
inline void swap(blocking_concurrent_queue<T, Traits>& a, blocking_concurrent_queue<T, Traits>& b) MOODYCAMEL_NOEXCEPT
{
	a.swap(b);
}

}	// end namespace mcl
