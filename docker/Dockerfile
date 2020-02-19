FROM python:3.7

ENV ARGS=''

RUN pip install nanome numpy scipy

COPY . /app
WORKDIR /app

CMD python run.py ${ARGS}
