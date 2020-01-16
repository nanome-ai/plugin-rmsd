FROM python:3.7

ENV PLUGIN_SERVER=plugins.nanome.ai

COPY . /app
WORKDIR /app

RUN pip install nanome
RUN pip install numpy
RUN pip install scipy

CMD python -m nanome_rmsd.RMSD -a ${PLUGIN_SERVER}