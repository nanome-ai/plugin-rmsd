FROM python:3.7

ENV ARGS=''

RUN pip install numpy scipy

ARG CACHEBUST
RUN pip install nanome

COPY . /app
WORKDIR /app

CMD python run.py ${ARGS}
