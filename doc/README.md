---
title: "Thesis Proposal"
author: "Muhammad Elsadany"
output:
  pdf_document: default
  word_document: default
---

### **1. Introduction**


*Background and Significance*

Genetics has been one of the most growing science fields, and yet we still know too little about diseases mechanisms or be able to provide proper treatments. Individual's genetics has been known to be a major player in differentiating response to different pharmacological compounds. The field of pharmacogenomics is rapidly growing, with the aim of providing personalized interventions for patients by uncovering the genetic basis of inter-individual variability in drug response. However, currently available pharmacogenomic services are limited in their ascertainment of genetic variation, with an emphasis on metabolism and known drug targets and mechanisms of action. This approach may result in blind spots relating to off-target effects as well as additional therapeutic mechanisms. A more comprehensive understanding of drug response could provide recommendations that maximize therapeutic effect while minimizing side effects.

In particular, there is a critical need for personalized drug recommendations for patients with psychiatric disorders. These patients often have complex medication regimens and can experience a wide range of responses to treatment, from total remission to life-threatening side effects. By integrating whole-genome genotype information, expression quantitative trait loci (eQTL), Polygenics Scores (PGS) of psychiatric traits, and transcriptional perturbation assays, it may be possible to develop models with reduced bias that can recommend drugs on an individual basis.

This research proposal aims to address the critical need for personalized drug recommendations for patients with psychiatric disorders by integrating disparate data sources to develop more comprehensive models of drug response. The proposed research builds on previous studies that have identified genetic and transcriptional signatures associated with psychiatric disorders, as well as studies that have investigated the effects of small molecules on gene expression. By integrating these data sources, we aim to develop a more comprehensive understanding of drug response that takes into account the complex genetic and environmental factors that contribute to inter-individual variability.

The significance of this research lies in its potential to improve patient outcomes by providing personalized drug recommendations that maximize therapeutic effect while minimizing side effects. By developing more comprehensive models of drug response that take into account the complex genetic and environmental factors that contribute to inter-individual variability, we hope to provide clinicians with a powerful tool for optimizing patient care. Ultimately, the proposed research could lead to improved treatment outcomes and a better quality of life for patients with psychiatric disorders.


*Research Question*

This study aims to address the following research questions:

* Can integrating whole-genome genotype information, expression quantitative trait loci (eQTL), and transcriptional perturbation assays improve the accuracy of personalized drug recommendations?

* Can the integration of high-resolution brain gene expression data and drug perturbation signatures help identify spatial patterns in susceptibility to small molecules in the brain, and can these patterns be used to further prioritize or de-prioritize drug recommendations based on anticipated effects in specific brain regions?

* Can the development of a new word embedding space capture the subjective experience of taking psychoactive medication, and can it be used to improve personalized drug recommendations?


*Objectives of the study* 

*this is redundant with the section above, probably only keep one?*

The objectives of this study are to:

* Integrate whole-genome genotype information, eQTL, GWAS, and transcriptional perturbation assays to recommend drugs on an individual basis that “normalize” disease-associated transcriptional signatures.

* Integrate high-resolution brain gene expression data and trait maps with drug perturbation signatures to identify spatial patterns in susceptibility to small molecules in the brain and use them to further prioritize or de-prioritize drug recommendations based on anticipated effects in specific brain regions.

* Develop a new word embedding space that captures the subjective experience of taking psychoactive medication and uses it to improve personalized drug recommendations.


*Hypothesis*

This study is based on the following hypotheses:

* Integration of whole-genome genotype information and eQTL data is sufficient enough to predict gene expression in different brain tissue, which could be integrated with transcriptional perturbation assays to predict personalized drug response. 

* Integration of high-resolution brain gene expression data with drug perturbation signatures can efficiently provide a gradient of drug activity in the brain. The drug gradient is hypothesized to be equivalent to anticipated effects in specific brain regions, that could be related to drug side effects.

* Psychoactive medications are more likely to cause different experiences that could not be captured by objective measures. We are hypothesizing that there is a highly similar subjective experiences among patients taking the same drug, and could not be measured through the first two approaches. 


### **2. Literature Review**

*overview of pharmacogenomics and personalized medicine*

Pharmacogenomics, also known as pharmacogenetics, is the branch of science that looks at how a person's genes influence how they react to pharmaceuticals. Its long-term objective is to assist physicians in choosing the medications and dosages that are ideal for every patient. It falls under the category of precision medicine, which tries to treat every patient uniquely. 

Pharmacogenomics, a rapidly developing field in personalized medicine, aims to elucidate how variations in genes can affect a patient's response to drugs. The influence of genes on drug metabolism and efficacy is well-established, as they encode for enzymes and proteins responsible for the breakdown and uptake of medications in the body.

Of particular interest are genes that encode for enzymes involved in drug metabolism, such as CYP2D6, which acts on a quarter of all prescription drugs. Multiple variations of this gene exist, with some individuals having multiple copies of it. These genetic variations can result in differences in enzyme activity, with some variants leading to a hyperactive enzyme that metabolizes drugs at a faster rate than normal. This can result in drug overdose, particularly in the case of codeine, which is metabolized by CYP2D6 to produce its active form, morphine. Conversely, some variants of CYP2D6 produce an enzyme that is non-functional or less active, leading to reduced or absent drug efficacy.

Therefore, understanding the impact of genetic variations on drug response is crucial in ensuring safe and effective drug therapy. Pharmacogenomics provides a valuable tool for predicting drug response based on an individual's genetic makeup, enabling personalized medicine approaches for improved patient outcomes.


*prev studies on drug response and pharmacogenomics*


* *points to be mentioned*

•	Mental disorders are a major challenge for individuals and society. Conditions such as schizophrenia, major depression, and anxiety disorders require long-term treatment with psychoactive drugs. Although there have been more than two-hundred drugs developed in the last six decades, they still can have variable effects between patients due to differences in drug metabolism and action. Consequently, increasing the dosage of medication does not necessarily lead to better treatment outcomes. 

•	Therapeutic drug monitoring (TDM) has emerged as a promising solution to these challenges, particularly for mood stabilizers, antidepressants, and antipsychotics. TDM has the potential to reduce variability, speed up clinical improvement, and improve drug tolerability and safety.

•	one of the main advantages of the proposed approach is that we are focusing on the entire genome, not only the genes involved in a drug response or metabolism.

•	Early treatment outcomes are frequently poor, with 30–50% of patients failing first-line antidepressant medication due to inefficiency or intolerance, according to estimates [2]. 

•	Moreover, about 25,000 people each year in the United States visit emergency rooms as a result of side effects brought on by antidepressants [3].

•	Before identifying a medication that reduces depressive symptoms with few side effects, patients frequently try with a variety of antidepressant regimens. The psychological and societal costs of repeatedly taking drugs that "do not work" can be distressing for the individual and highlight the need for better drug selection and dosage tactics because antidepressant pharmacotherapy studies frequently take a minimum of 6-8 weeks.

•	The FDA has released comments and warning letters on pharmacogenetic testing in response to concerns over the marketing of these tests. The efficacy of clinical pharmacogenetic testing may not be fully supported by clinical data, according to a safety communication released on November 01, 2018. The statement in this safety communication that "the relationship between DNA variations and the effectiveness of antidepressant medication has never been established" particularly emphasized the use of pharmacogenetic testing to guide antidepressant drug prescribing [4].

•	In keeping with the FDA goal to safeguard and advance the public's health, it's critical to act right away to make sure that the claims being made about the pharmacogenetic tests currently available are supported by reliable research. That can be done by taking measures that safeguard patients while also advancing the creation of analytically and clinically validated pharmacogenetic tests. Recently, the FDA released a new web-based resource that includes a table including some of the pharmacogenetic associations with a last update on October 26, 2022 [5]. Some of these have detailed information regarding therapeutic management, but the majority of the associations listed have not been assessed in terms of the effect of genetic testing on clinical outcomes, such as improved therapeutic effectiveness or increased risk of particular adverse events. This version of the table is restricted to pharmacogenetic associations linked to drug transporter, drug metabolizing enzyme, and gene variations associated with a susceptibility for certain adverse outcomes.

•	pharmacogenomic testing could be done in two forms: single-gene or a multi-panel testing. 

•	Most of the developed pharmacogenetic tests investigate the main drug gene targets, and the advanced ones with a multigene panel look into different genes that include the ones involved in the process of metabolizing the drug. Yet, none of these tests consider a whole-genome approach, which is one of our main focuses in the proposed study. 

•	The FDA now advises against using direct-to-consumer testing for making medical choices, even though there has been considerable debate about the use of pharmacogenomics in clinical practice. However, the FDA allowed the marketing of the 23andMe Personal Genome Service Pharmacogenetic Reports test as a direct-to-consumer test with special controls for informing discussions with a healthcare professional about genetic variants that may be related to a patient's capacity to metabolize some medications [6].



*gene expression, eQTL, and GWAS related to psychiatric traits*


*transcriptional perturbation assays and their application to drug discovery*



### **3. Methods**

*study design and sample selection*


*data sources and management*


*statistical analyses and models*

*ethical considerations*



### **4. Specific Aims**

*Aim 1*


*Aim 2*


*Aim 3*



### **5. Expected outcomes and Significance**

*Potential impact of the study on personalized medicine and drug discovery*


*Contribution to the understanding of the genetic basis of drug response and psychiatric traits*


*Implications for clinical practice and patient outcomes*



### **6. Timeline**

*estimated timeline for data analysis and publication*


### **7. Conclusion**

*summary of the proposed research and its significance*

*implications for future research and clinical practice*