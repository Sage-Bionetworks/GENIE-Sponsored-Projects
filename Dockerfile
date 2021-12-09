FROM python:3.8 

WORKDIR /usr/src/GENIE-Sponsored-Projects

COPY requirements.txt ./ 
RUN pip install --no-cache-dir -r requirements.txt 

RUN apt-get -y install git
RUN git clone git://github.com/cBioPortal/cbioportal.git ../cbioportal

COPY . .

ENTRYPOINT [ "python", "-m", "geniesp" ]
