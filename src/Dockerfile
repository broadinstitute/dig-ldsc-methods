FROM ubuntu:20.04

RUN apt-get update && \
    apt-get install -y git build-essential libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev curl wget zlib1g-dev python3-pip

WORKDIR /usr/src

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY . .


