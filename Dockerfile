FROM python:3.8 

WORKDIR /usr/src/GENIE-Sponsored-Projects

RUN git clone git://github.com/cBioPortal/cbioportal.git ../cbioportal

COPY . .
RUN pip install --no-cache-dir -r requirements.txt 

ENTRYPOINT [ "python", "-m", "geniesp" ]
