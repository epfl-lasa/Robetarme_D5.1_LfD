test_lagsDS can run

Steps 1-5 are to fit the GMM and learn the global DS. Steps 6-7 corrrespond to the lags-learning part.

Step 6: You can see the Gaussians that are fitted on each linear region
Then in this variable
choosen_active = 1:K;
 is where you define which regions should be locally activated; i.e. which Gaussians.

Line 380 sets the width of the RBF kernel in the activation function, which defines the width/strength of the locally active regions.

Step 7: Here is where I estimate the parameters of the locally active regions. Option 1 and 2 are two different implementations of the optimization problem in my thesis. They will both give you stable, yet conservative, estimates.  With Option 3 you can manually define the tracking factor of the local systems (i.e. eigenvalue ratio). This might be useful for you as you would like to define this ratio from the impedance learned from the human.

Step 8: To visualise the 2D LAGS-DS