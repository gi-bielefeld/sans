FROM alpine:latest
RUN apk --no-cache add build-base zlib-dev
WORKDIR /app
ENV PATH="$PATH:/app"
COPY . .  
RUN make

CMD ["/bin/sh",  "-c", "/app/SANS"]
