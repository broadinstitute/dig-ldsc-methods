FROM ubuntu:20.04

RUN apt-get update && apt-get install -y unzip python3-pip

WORKDIR /usr/src

COPY requirements.txt .
RUN pip3 install --no-cache-dir -r requirements.txt
RUN pip3 --no-cache-dir install --upgrade awscli

COPY src .
