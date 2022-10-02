# Image Classification with Pseudo Labels 

### Visualize Dataset 
This dataset consists of 11 categories of food. For this scenario, we only have 
~30% of our dataset labelled. Let's investigate how pretraining and self-training 
can help with our situation.
1. Training set (with labels):  5% (Specifically, 493 images, 
                                  i.e. on average ~40 for each category)
2. Training set (unlabeled)  : 70%
3. Validation set            : 25%

![input](imgs/input.png)

### Without Pretraining
![input](imgs/baseline_40.png)

### With Pretraining
![input](imgs/pretrain_40.png)

**Pseudo Labels**  
![Alt Text](pseudo.png)
