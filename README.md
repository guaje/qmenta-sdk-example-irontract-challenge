# QMENTA SDK Tool for the IronTract Challenge Round II

Repository cloned from [QMENTA SDK Tool Example](https://github.com/qmentasdk/qmenta-sdk-example).

## Contents

### Dockerfile

This file contains the sequence of instructions to build a new tool image. You can setup environment paths, run commands during the image building stage and copy files (see [Dockerfile commands](https://docs.docker.com/get-started/part2/)).

For this particular application [DIPY](https://dipy.org/) and some of its [dependencies](https://dipy.org/documentation/1.2.0./dependencies/#dependencies) are installed.

### Tool script

The main script that is executed when a tool is launched on the [QMENTA platform](https://platform.qmenta.com/). This script typically performs the actions shown below using the [QMENTA SDK](https://docs.qmenta.com/sdk) functions where suitable:

1. Download the input data from the [IronTract challenge](https://irontract.mgh.harvard.edu/) to the container.
2. Process it using [DIPY](http://dipy.org).
3. Upload the results.

## Build the tool image

Use [Docker](https://www.docker.com/get-docker) to build a new image using the Dockerfile:
~~~~
docker build -t image_name .
~~~~
Where `image_name` should conform to the syntax `my_username/my_tool:version`.

> The first build may take several minutes since it will need to generate the image layer containing the software installation.

Alternatively, take a look at the `standalone.Dockerfile` to see how to install the SDK in an image based on Ubuntu.

### Report template

This tool code uses this HTML template to populate some fields with the patient data and the analysis results to generate a PDF report.

## Test the tool locally

Optionally, the `test_tool.py` script can be used to locally launch your tool image if you specify the input files and the required values for you settings (see `settings_values.json`):
~~~~
mkdir analysis_output

python test_tool.py image_name example_data analysis_output \
    --settings settings.json \
    --values mock_settings_values.json
~~~~

## Add the tool to the [QMENTA platform](https://platform.qmenta.com/)

To add a tool image to your list of tools you will first need to login push it to your [Docker Hub](https://hub.docker.com/) registry:
~~~~
docker login

docker push image_name
~~~~
To register the tool to the [QMENTA platform](https://platform.qmenta.com/), access the Analysis menu, and go to My Tools. You will need to enter your credentials of your registry, the name and a version number for the tool, its settings configuration file and a description. You can find an example of the settings configuration file in this repository (`settings.json`).
