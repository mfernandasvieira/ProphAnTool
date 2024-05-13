# ProphAnTool (Prophage Analysis Tool)
ProphAnTool is a comprehensive tool designed for the efficient search and analysis of large bacterial genome datasets. 
It automates the process, allowing users to swiftly extract insights from their genomic data. 
This tool integrates various publicly available tools into a unified pipeline, streamlining the analysis process and providing all results in one consolidated output.

## Installation

To install ProphAnTool, follow these steps:

1. Clone Repository: Clone the ProphAnTool repository to your local machine.

```
git clone https://github.com/mfernandasvieira/ProphAnTool.git

```

2. Set Up Environment: Navigate to the cloned repository and create the Conda environment using the provided environment configuration file.

```
conda env create -f environment.yml
```

3. Activate Environment: Activate the Conda environment.

```
conda activate prophantool
```

## Usage
Once the environment is set up, you can start using ProphAnTool. Follow these instructions:

1. Navigate to ProphAnTool Directory: Move into the ProphAnTool directory.
```
cd ProphAnTool
```

2. Run ProphAnTool: Execute ProphAnTool by running the following command:
```
python prophantool.py -s Snakefile -c default_config.json
```

This command initiates the ProphAnTool pipeline using the specified Snakefile and configuration file. Make sure to adjust the filenames if necessary.

3. Access Results: After the analysis is complete, you can access the generated results in the output directory called results.

## License
ProphAnTool is licensed under the MIT License. See the LICENSE.md file for more details.