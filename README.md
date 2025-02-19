# Core Switching Circuit
This is the supporting code to the paper titled **NetCF : A network control-based framework to reveal the molecular mechanism of phenotype switching in lung cancer cells**

The Core Switching Circuit algorithm is a part of the NetCF. It identifies the essential circuit that explains the effects when a specific drug and its combination target are regulated simultaneously. This facilitates the analysis of the molecular mechanism by which the drug and its combination target act on cells. For the front part of the NetCF framework, please refer to [NetCF](https://github.com/namheee/NetCF).  


## Installation
No special installation is required. Simply clone this repository to a suitable location and run `main.py` using Python 3 as follows: 

```bash
git clone https://github.com/yena2bell/CoreSwitchingCircuit
cd CoreSwitchingCircuit
python3 main.py command_file_address
```

For required argument (command_file_address), please refer to the Basic Usage Example section.

### Requirements
NumPy (v1.25.0) https://numpy.org/

Matplotlib (v3.7.1) https://matplotlib.org/

## Basic usage example
To run this algorithm, the following data must be prepared:  

1. A TSV file containing network structure information. The format can be referenced in `example_data/links_of_KRAS_mutant_lung_cancer_cells_model.tsv`.  
2. Information on the nodes affected by the drug within the given network structure.  
3. Information on the nodes affected by the combination target within the given network structure.  
4. Information on the nodes that determine the phenotype within the given network structure.  
5. A TSV file containing the average node activity profile under nominal conditions.  
6. A TSV file containing the average node activity profile when only the drug is applied.  
7. A TSV file containing the average node activity profile when only the combination target is applied.  
8. A TSV file containing the average node activity profile when both the drug and combination target are applied.  

Data for points 5–8 can be computed using **Notebook 4** from [NetCF](https://github.com/namheee/NetCF).  
The required format can be referenced in `example_data/average_node_activity_profile_for_KRAS_ON.tsv`.  

Once all required data is prepared, the information should be recorded in a **command file**.  
An example command file is provided in `command_file_example.txt`.  
The `save address` field in the command file specifies the output file where the results will be stored.  

After setting up the command file, run `main.py` with the first argument as the command file address, as shown below:  

```bash
python3 main.py command_file_address
```

### Required Parameters  

The following parameters are required to run this algorithm:  

- **Threshold for the initially perturbed region of the drug target**  
- **Threshold for the initially perturbed region of the combination target**  
- **Threshold for the feedback score**  
- **Threshold for significant regulator determination**  
- **Maximum feedback length to search**  

When executing `main.py` with Python, a command-line interface (CLI) prompts the user to enter each parameter.  
During this process, **matplotlib-based figures** will **pop up**, providing useful references for parameter selection.  

For details on how to determine each parameter and utilize the displayed figures, please refer to **Supplementary XX** of the paper.  


### Output Format
An example of the output is as follows: 

```
initially perturbed region
{'JUN', 'MYC', 'SMAD2', 'SP1', 'SMAD3', 'MAPK3', 'CREB1', 'TP53', 'MAPK14'}

0th step
PPR:{'JUN', 'MYC', 'SMAD2', 'SP1', 'SMAD3', 'CEBPA', 'MAPK3', 'CREB1', 'TP53', 'MAPK14'}
added feedback:(('MYC', '-', 'CEBPA'), ('CEBPA', '-', 'MYC'))
Ftb of the feedback:0.33333333333333337
feedback score of the feedback:0.7741787823402693

1th step
PPR:{'JUN', 'MYC', 'SMAD2', 'E2F1', 'SP1', 'SMAD3', 'CEBPA', 'MAPK3', 'CREB1', 'TP53', 'MAPK14'}
added feedback:(('MYC', '+', 'SP1'), ('SP1', '+', 'E2F1'), ('E2F1', '+', 'MYC'))
Ftb of the feedback:0.7142857142857143
feedback score of the feedback:0.7233366816659949

...

10th step
PPR:{'MYC', 'STK11', 'CDK1', 'STAT3', 'ESR1', 'CDK2', 'CDKN1A', 'CREB1', 'PRKAA1', 'MAPK14', 'JUN', 'MTOR', 'SMAD2', 'PPP1CA', 'E2F1', 'SP1', 'MAPK1', 'SMAD3', 'RPS6KB1', 'MAPK3', 'GSK3B', 'CEBPA', 'TP53', 'AKT3'}
added feedback:(('PRKAA1', '-', 'MTOR'), ('MTOR', '+', 'RPS6KB1'), ('RPS6KB1', '-', 'GSK3B'), ('GSK3B', '-', 'MYC'), ('MYC', '-', 'CDKN1A'), ('CDKN1A', '-', 'CDK1'), ('CDK1', '-', 'PPP1CA'), ('PPP1CA', '-', 'MAPK1'), ('MAPK1', '+', 'MAPK3'), ('MAPK3', '-', 'STK11'), ('STK11', '+', 'PRKAA1'))
Ftb of the feedback:0.5
feedback score of the feedback:0.6584334002226464
```

Here, the **initially perturbed region (IPR)** represents the union of the IPRs of the drug and the combination target.  

During the identification of the **core switching circuit**, the **Perturbation-Propagated Region (PPR)** is iteratively expanded. In the output above, the **i-th step** corresponds to the information when the PPR has been expanded for the **i-th time**.  

- **PPR**: The expanded **Perturbation-Propagated Region** at the **i-th step**.  
- **Ftb of the feedback**: The **Feedback Transition-Barrier (Ftb)** value calculated for the feedback added at the **i-th step**. A **higher Ftb value** indicates that the activation of this feedback is more difficult.  
- **feedback score of the feedback**: The **feedback score** calculated for the feedback added at the **i-th step**. A **higher feedback score** indicates that nodes significantly affected by the drug and combination target are included in this feedback.  

