FROM python:3.8 

WORKDIR /usr/src/GENIE-Sponsored-Projects

COPY . .
RUN pip install ./

RUN git clone https://github.com/cBioPortal/cbioportal.git -b v5.3.19 ../cbioportal
