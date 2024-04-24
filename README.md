# Panama_Antibiotics

![image info](Results/overall_composition.png)

The vast majority (620) of all retained ASVs show no difference between the initial antibiotic dosed and undosed treatments (`fdr_preDose`) nor among any of the post disease exposure treatment combinations (`fdr_treatment`) including the healthy outcome antibiotic treated, healthy outcome not antibiotic treated, and disease outcome not antibiotic treatment. We used three planned post-hoc tests to identify ASVs of interest. Specifically, is there more of the ASV in the disease outcome than the average of the two healthy outcomes regardless of the antibiotic treatment (`fdr_disease.v.avg`). The second and third tests individually test that the no antibiotic treated disease outcome has more of the ASV than the no antibiotic healthy outcome (`fdr_disease.v.h`) and than the antibiotic treated healthy outcome (`fdr_disease.v.anti`).

Of the ASVs which do show some response 48 show differences between treatments but do not show any of the specific *a priori* differences expected. A total of 6 show significant differences between treatments and all three *a priori* contrasts with a further 2 also being significantly influenced by the pre-dose antibiotic treatment.
![image info](Results/asv_upset.png)

These 8 samples include *Cysteiniphilum litorale* along with five Roseobacteraceae, one Rhodobacterales, one Pseudomonadota, and *Marisediminitalea aggregata*. One of the Roseobacteraceaes, the Rhodobacteralesm and the *M. aggregata* all were significantly affected by the initial antibiotic treatment prior to disease homogenate dosing.
![image info](Results/asvs_changing_postExposure.png)

Complete ASV results are [here](Results/individual_asv_results.csv)
