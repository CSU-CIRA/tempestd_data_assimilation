FROM pyncepbufr as base

MAINTAINER Jim Fluke <james.fluke@colostate.edu>

#Health check port 
EXPOSE 5000

WORKDIR /app

COPY . /app

ENTRYPOINT [ "python", "multiproc_converter.py" ]
