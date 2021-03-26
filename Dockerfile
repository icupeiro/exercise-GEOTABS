FROM michaelwetter/ubuntu-1604_jmodelica_trunk:latest

RUN pip install --user --no-cache --upgrade pip && \
RUN pip install --user --no-cache notebook

ENV ROOT_DIR /usr/local
ENV JMODELICA_HOME $ROOT_DIR/JModelica
ENV IPOPT_HOME $ROOT_DIR/Ipopt-3.12.4
ENV SUNDIALS_HOME $JMODELICA_HOME/ThirdParty/Sundials
ENV SEPARATE_PROCESS_JVM /usr/lib/jvm/java-8-openjdk-amd64/
ENV JAVA_HOME /usr/lib/jvm/java-8-openjdk-amd64/
ENV PYTHONPATH $PYTHONPATH:$JMODELICA_HOME/Python:$JMODELICA_HOME/Python/pymodelica
ENV PATH="/home/developer/.local/bin:${PATH}"

USER root

ARG NB_USER
ARG NB_UID
ENV USER ${NB_USER}
ENV HOME /home/${NB_USER}

RUN adduser --disabled-password \
    --gecos "Default user" \
    --uid ${NB_UID} \
    ${NB_USER}

WORKDIR $HOME
USER ${USER}

RUN apt-get update && \
	apt-get install -y git



RUN pip install --user future
RUN pip install --user pandas
RUN pip install --user scipy
RUN pip install --user numpy
RUN pip install --user matplotlib
RUN pip install --user ipykernel==4.7.0
RUN pip install --user pygfunction

COPY Introduction.ipynb $HOME
COPY PartI_TABS.ipynb $HOME
COPY PartII_StaticCalc.ipynb $HOME
COPY PartIII_EnergySystemDesign.ipynb $HOME
COPY PartIV_DynamicSim.ipynb $HOME
COPY PartV_hybridGEOTABS.ipynb $HOME
COPY loadCalc.xlsx $HOME
COPY resistanceCalculator.py $HOME
COPY fig $HOME/fig
COPY results $HOME/results
