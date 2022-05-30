FROM continuumio/miniconda3
WORKDIR /IntegrationSiteSeq
COPY requirements.txt .
RUN conda create --name integsite --file requirements.txt
RUN conda activate integsite
COPY IntegrationSite.py .
ENTRYPOINT ["python", "integsite"]