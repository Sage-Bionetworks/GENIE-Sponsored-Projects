FROM python:3.8 

WORKDIR /usr/src/app 

COPY .synapseConfig ./
COPY cbioportal/ ./
COPY requirements.txt ./ 
RUN pip install --no-cache-dir -r requirements.txt 

COPY . .

ENTRYPOINT [ "python", "-m", "geniesp" ]
