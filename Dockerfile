FROM python:3.9-slim
RUN apt-get update && apt-get install -y build-essential zlib1g-dev bash default-jre xvfb && apt-get clean
RUN pip install --no-cache-dir pygtrie cairosvg
WORKDIR /app
ENV PATH="$PATH:/app:/app/scripts"

ADD --checksum=sha256:8680653437a711750dfbc9a1a812b6069b0570badb8b78ee903c0e2cdfb8cc75 --chmod=777 https://software-ab.cs.uni-tuebingen.de/download/splitstree4/splitstree4_unix_4_19_2.sh /tmp/splitstree.sh
RUN { echo 'o' ; echo '1' ; echo '/opt/splitstree' ; echo 'X,2,3' ; echo "n" ; echo '5' ; echo '2000' ; } | /tmp/splitstree.sh
ENV DISPLAY=":0.0"
RUN ln -s /opt/splitstree/SplitsTree /usr/bin/SplitsTree
COPY . .
RUN make
RUN chmod 777 -R /app
