Looking into emu results, we are seeing cases where the same read can match to reference sequences with two different taxids.

sample	| 473972_uncultured_Metastrongylid | 1258553_Varestrongylus_eleguneniensis
--------|----------------------------------|--------------------------------------
sample3 | 846.680954859078 | 44.3190451409212

```grep "473972" pid_debug_output.txt | head```
d1c4d768-14df-4dcb-b539-ce0c0cca9a49	473972:nem_its2:1697	1	99.80  
d1c4d768-14df-4dcb-b539-ce0c0cca9a49	473972:nem_its2:1693	2	99.60  
d1c4d768-14df-4dcb-b539-ce0c0cca9a49	473972:nem_its2:1683	3	99.39  
d1c4d768-14df-4dcb-b539-ce0c0cca9a49	473972:nem_its2:1695	3	99.39  
d1c4d768-14df-4dcb-b539-ce0c0cca9a49	473972:nem_its2:1692	15	96.96  
d1c4d768-14df-4dcb-b539-ce0c0cca9a49	473972:nem_its2:1691	10	97.98  
d1c4d768-14df-4dcb-b539-ce0c0cca9a49	473972:nem_its2:1684	3	99.35  
9d941bd5-2b60-4caa-97f3-fd77794d217e	473972:nem_its2:1697	1	99.80  
9d941bd5-2b60-4caa-97f3-fd77794d217e	473972:nem_its2:1693	2	99.60  
9d941bd5-2b60-4caa-97f3-fd77794d217e	473972:nem_its2:1683	3	99.39  

```grep "1258553" pid_debug_output.txt | head ```
d1c4d768-14df-4dcb-b539-ce0c0cca9a49	1258553:nem_its2:1694	2	99.60  
d1c4d768-14df-4dcb-b539-ce0c0cca9a49	1258553:nem_its2:5503	3	99.39  
d1c4d768-14df-4dcb-b539-ce0c0cca9a49	1258553:nem_its2:5505	3	99.39  
d1c4d768-14df-4dcb-b539-ce0c0cca9a49	1258553:nem_its2:1696	3	99.39  
d1c4d768-14df-4dcb-b539-ce0c0cca9a49	1258553:nem_its2:5504	5	98.99  
d1c4d768-14df-4dcb-b539-ce0c0cca9a49	1258553:nem_its2:10728	2	99.55  
9d941bd5-2b60-4caa-97f3-fd77794d217e	1258553:nem_its2:1694	2	99.60  
9d941bd5-2b60-4caa-97f3-fd77794d217e	1258553:nem_its2:5503	3	99.39  
9d941bd5-2b60-4caa-97f3-fd77794d217e	1258553:nem_its2:1696	3	99.39  
9d941bd5-2b60-4caa-97f3-fd77794d217e	1258553:nem_its2:5505	3	99.39  

The debug output shows that the reads can be matched to either reference sequence with similar percentage identity.
We wanted to be able to see percentage ID of the reads and be able to have more flexibility with the output, including getting haplpotype sequences that we could examine and carry out donstream analysis of.
So even though emu is really fast computationally, we have decided to switch to a haplotype-based approach, where we generate and polish haplotypes, and blast the polished haplotypes against databases.
This process is similar to how we generated Setaria and Onchocerca reference sequences, and references codes from Gandasegui (https://github.com/Gandasegui/Incognita_reanalysis/blob/main/01_scripts/04_generating_haplotype_counts.md).
