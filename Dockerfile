FROM python:3.10-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    STREAMLIT_SERVER_ADDRESS=0.0.0.0 \
    STREAMLIT_SERVER_PORT=8501

WORKDIR /app

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        build-essential \
        gcc \
        g++ \
        git \
        libxml2-dev \
        libsbml5-dev \
    && rm -rf /var/lib/apt/lists/*

COPY pyproject.toml README.md requirements.txt ./
COPY .streamlit ./.streamlit
COPY src ./src
COPY app ./app
COPY pure_most_detailed.py fMGG_synthesis_parameters_CRN.csv PURE_TXTL_initial_values_Final.csv ./

RUN python -m pip install --upgrade pip setuptools wheel \
    && python -m pip install --no-cache-dir "numpy==1.26.4" "Cython<3" \
    && python -m pip install --no-cache-dir --no-build-isolation "bioscrape>=1.2" \
    && python -m pip install --no-cache-dir -r requirements.txt

EXPOSE 8501

CMD ["streamlit", "run", "app/streamlit_app.py"]
