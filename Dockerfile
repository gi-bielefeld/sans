FROM ubuntu:22.04
RUN apt-get update && apt-get install -y build-essential zlib1g-dev && apt-get clean
WORKDIR /app
ENV PATH="$PATH:/app"
COPY . .  
RUN make

CMD ["/bin/sh",  "-c", "/app/SANS"]
