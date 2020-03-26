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

USER developer

RUN git clone https://github.com/open-ideas/IDEAS.git $HOME/IDEAS
RUN git clone https://github.com/icupeiro/model-GEOTABS.git
RUN pip install --user --no-cache-dir notebook==5.*
RUN pip install --user pandas
RUN pip install --user flask
RUN pip install --user ipykernel==4.7.0
RUN pip install --user pygfunction

COPY JModelica_test.ipynb $HOME
