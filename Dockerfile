FROM michaelwetter/ubuntu-1604_jmodelica_trunk

ENV ROOT_DIR /usr/local
ENV JMODELICA_HOME $ROOT_DIR/JModelica
ENV IPOPT_HOME $ROOT_DIR/Ipopt-3.12.4
ENV SUNDIALS_HOME $JMODELICA_HOME/ThirdParty/Sundials
ENV SEPARATE_PROCESS_JVM /usr/lib/jvm/java-8-openjdk-amd64/
ENV JAVA_HOME /usr/lib/jvm/java-8-openjdk-amd64/
ENV PYTHONPATH $PYTHONPATH:$JMODELICA_HOME/Python:$JMODELICA_HOME/Python/pymodelica
ENV PATH="/home/developer/.local/bin:${PATH}"

USER root

RUN apt-get update && \
	apt-get install -y git

 
ENV HOME /home/developer
WORKDIR $HOME

ARG NB_USER=jovyan
ARG NB_UID=1000
ENV USER ${NB_USER}
ENV NB_UID ${NB_UID}
ENV HOME /home/${NB_USER}

USER developer
--uid 1000

RUN pip install --user --no-cache-dir notebook==5.*
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
