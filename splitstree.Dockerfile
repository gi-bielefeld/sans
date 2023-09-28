FROM python:3.9-slim
RUN apt-get update && apt-get install -y wget bash default-jre xvfb
RUN pip install --no-cache-dir pygtrie

RUN wget -nv -O /tmp/splitstree.sh https://software-ab.cs.uni-tuebingen.de/download/splitstree4/splitstree4_unix_4_19_2.sh
RUN chmod 777 /tmp/splitstree.sh
RUN { echo 'o' ; echo '1' ; echo '/opt/splitstree' ; echo 'X,2,3' ; echo "n" ; echo '5' ; echo '2000' ; } | /tmp/splitstree.sh
ENV DISPLAY=":0.0"
RUN ln -s /opt/splitstree/SplitsTree /usr/bin/SplitsTree
COPY ./scripts/ /usr/bin
