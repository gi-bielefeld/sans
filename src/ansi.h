#define $out  std::cout
#define $log  std::cerr

#define $err  std::cerr << "\33[38;5;131m"
#define $warn std::cerr << "\33[38;5;137m"
#define $note std::cerr << "\33[38;5;132m"

#define $link std::cerr << "\33[38;5;67m"
#define $lite std::cerr << "\33[38;5;244m"
#define $_              << "\33[2K\r"

#define     $              std::flush
#define    _$  "\33[0m" << std::flush
#define  end$              std::endl
#define _end$  "\33[0m" << std::endl
