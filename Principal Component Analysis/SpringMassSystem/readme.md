This folder contains the analysis of the dynamics of spring-mass system using principal component analysis (PCA) and singular value decomposition (SVD).

Different caees were taken here:

1. Case 1: this is the ideal case where you consider a small displacement of the mass in the z direction and the ensuing oscillations. In this case, the entire motion is in the z direction with the simple harmonic motion being observed (camN_1.mat where N=1,2,3).

2. Case 2: this is the noisy case where you repeat the ideal case experiment, but this time, introduce camera shake into the video recording. This should make it more difficult to extract the simple harmonic motion. But if the shake isn't too bad, the dynamics will still be extracted with the PCA algorithms. (camN_2.mat where N=1,2,3)

3. Case 3: in this case, you have a horizontal displacement. Here, the mass is released off-center to produce motion in the x-y plane as well as the z-direction. Thus there is both a pendulum motion and simple harmonic oscillations. See what the PCA tells us about the system. (camN_3.mat where N=1,2,3)

4. Case 4: in this case, you have a horizontal displacement and rotation. Here, the mass is released off-center and rotates to produce motion in the x-y plane, rotation as well as the z-direction. Thus there is both a pendulum motion and simple harmonic oscillations. See what the PCA tells us about the system. (camN_4.mat where N=1,2,3)

The datasets of all these different cases were taken from: [Link](https://drive.google.com/drive/folders/1SQ77P5t5RUWCSucmk4jPFbufFMX8VrJG)
