# Re-emergence of respiratory syncytial virus following the COVID-19 pandemic in the United States: a modeling study
This repo contains the codes for https://www.medrxiv.org/content/10.1101/2021.07.19.21260817v1
1) Download parameters and codes in SimpleModel file. Run the simple-age-structure model to generate simulated epidemic data. 
2) If you have inpatient data or birth cohort observations, using the last chunk of code in simple-age-structure model to fit to your observation.
3) Using the parameters you get and run the forward simulation by changing the code 
results <- ode(y=yinit.vector, t=my_times,  
               func=simple_model, # your new models. The three .R files in the main page corresponding to three scenarios in the paper
            parms=parms ) # your new parameter set
