
1- ***Stability of Time-Stepping Schemes*** : Basic time-stepping schemes, e.g., forward and backward Euler, involve discretisation in both space and time for time-dependent partial differential equations. Issues of numerical stability as the solution is propagated in time and accuracy from the space-time discretisation are of great concern for any implementation. The optimal choice of the time-stepping scheme depends on its purpose to obtain a time-accurate discretisation of a highly dynamic flow problem or to march the numerical solution to a steady state starting with some reasonable initial guess. Every scheme eventually leads to an iteration procedure which the computer can use to advance the solution in time. This project assesses the credibility of the solution and determines the ratio (or a number) that controls the accuracy and stability of the chosen scheme. Equations such as wave and heat differential equations are to be examined using MATLAB software.

-------------------------------------------------------------

2- ***Independent Component Analysis for Image Separation***: Independent component analysis (ICA) is a statistical and computational technique for revealing hidden factors that underlie sets of random variables, measurements or signals. ICA can be seen as an extension to principal component analysis (PCA). What distinguishes ICA from other methods is that it looks for components that are both statistically independent, and non-Gaussian. Its applications are wide and varied and the data analysed by ICA could originate from many different application fields, including digital images and document databases, as well as economic indicators and psychometric measurements. The primary focus of this project will be the application of ICA in the context of digital image separation. The method proposed here is to separate two images relies on reversing the action of the singular value decomposition (SVD) method on two statistically independent images using MATLAB software.

-------------------------------------------------------------

3- ***Principal Component Analysis for face recognition***: Dimensionality reduction, typically achieved through a principal component analysis (PCA) or orthogonal mode decomposition (POD), has achieved remarkable success in providing a mathematical framework which is much more amenable to analysis. This has allowed a better characterisation of the physics, engineering or biological system of interest. PCA is an unsupervised learning method, which seeks to find natural clustering of data in a lower dimensional feature space, whose basis contains linearly, ordered, independent, orthonormal and uncorrelated vectors that correspond to the maximum-variance directions in the original space. Interestingly, PCA can be used in the context of the computer vision problem of human face recognition. The approach of using eigenfaces for recognition functions by projecting face images onto a feature space that spans the significant variations among known face images. The significant features are known as “eigenfaces”, because they are the eigenvectors (principal components) of the set of faces. The projection operation characterises an individual face by a weighted sum of the eigenface features, and so to recognise a particular face it is necessary only to compare these weights to those of known individuals. This project uses MATLAB software to examine this method.

-------------------------------------------------------------

4- ***Numerical Analysis of Shallow-Water Equations (Solution to the two-dimensional incompressible Navier-Stokes Equations with vorticity-streamfunction Approach.)***: Numerical simulations of the Navier-Stokes equations for studying two-dimensional flows of incompressible viscous fluid are generally based upon a primitive variables formulation velocity and pressure or vorticity-streamfunction formulation. The major difficulty arising with the former formulation comes from the coupling of the pressure with the velocity, to satisfy the incompressibility condition. The continuity equation contains only velocity components, and there is no direct link with the pressure as it happens for compressible flow through the density (the lack of evolution equation for the pressure in primitive variables formulation is the source of difficulty). The use of a vorticity-streamfunction formulation of the equations avoids this problem. The problem can be further simplified using the limit of Shallow Water Equations (SWE) where the horizontal length scale is much greater than the vertical length scale. SWE are of fundamental interest in several contexts. In one sense, the ocean can be thought of as a shallow-water description over the surface of the Earth. The atmosphere can be also thought of as a relatively thin layer of fluid (gas) above the surface of the Earth. To understand the motion in the shallow-water limit, in this project, we will investigate the concept of a stream-function ψ and vorticity ω formaulation as dependent variables and how they are coupled to one another in two dimensions and how to solve these equations numerically via MATLAB software using three methods:  
1-Finite Difference Method  
2-Fourier Transform  
3-Chebychev Spectral Differentiation   
A real-life application for these equations would be: "Finite Element Analysis of a Flow Over an Airfoil" (mentioned here as an individual project).

-------------------------------------------------------------

5- ***Numerical Analysis of Heat Equation***: The diffusion equation is a partial differential equation which describes density fluctuations in a material undergoing diffusion. This equation is also called the heat equation and also describes the distribution of
a heat in a given region over time. This equation can describes many natural science and engineering phenomena that involve diffusion, such as: the diffusion of heat in the skin as a result of a burning accident, the diffusion of electrons when their density in solids is not in equilibrium, the diffusion of the plasma based on the strength of an external magnetic field, and many more. This equation can be solved either analytically or numerically. It is mathematically relatively simple and analytical solutions can be easily found if geometry and boundary conditions are not too complicated. For the numericaly treatment, a plethora of methods can be used to solve such essential and basic equation. In this project, we will solve the 1D and 2D heat transfer equation, using three methods:   
1-Finite Difference Method  
2-Fourier Transform  
3-Chebychev Spectral Differentiation   

-------------------------------------------------------------

6- ***Image Denoising***: With the explosion in the number of digital images taken every day, the demand for more accurate and visually pleasing images is increasing. However, the images captured by modern cameras are inevitably degraded by noise, which leads to deteriorated visual image quality. Therefore, work is required to reduce noise without losing image features (edges, corners, and other sharp structures). Image denoising is to remove noise from a noisy image, so as to restore the true image. Image denoising plays an important role in a wide range of applications such as image restoration, visual tracking, image registration, image segmentation, and image classification, where obtaining the original image content is crucial for strong performance. While many algorithms have been proposed for the purpose of image denoising, in this project, we present two image denoising techniques, filtering and smoothing. For filtering, we will present two filers, Gaussian (bell shaped) and Shannon (square function) and decide which one is the better filter and why. However for smoothing, we will use the diffusion equation to smooth out an image. 

-------------------------------------------------------------

7- ***Time-Frequency Analysis and Gabor Transform***: Fourier transform analysis has long been recognised as the great tool for the study of analysing the properties of signals. Fourier Transform represents a signal in its frequency domain with no ability to view the distribution of various frequencies in time. This is one of the shortfalls of the Fourier transform analysis since it precludes a very practical notion, namely the notion of frequencies changing with time. One of the earliest mathematical techniques used to remedy this deficiency is the Gabor transform. It permits a time-frequency representation of a signal, hence making it possible to view the frequency spectrum locally in time. Indeed, this method is used extensively for analyzing speech and vocalisation patterns. For such applications, it is typical to produce a spectrogram that represents the signal in both the time and frequency domains. This project focuses on understanding this method and how it is applied on different types of signals using MATLAB software.

-------------------------------------------------------------

8- ***Classical and Quantum Oscillators***: Oscillators are physical systems whose evolution over time varies repeatedly around a central state of equilibrium. Oscillating cycles that are more or less periodic are found in all sectors of science, from quantum physics to cell biology, sociology or cosmology. In this project, we will look into studying the behaviour of three types of oscillators:   
- *Damped Harmonic Oscillator:* In classical mechanics, a harmonic oscillator is a system that, when displaced from its equilibrium position, experiences a restoring force F proportional to the displacement x. If F is the only force acting on the system, the system is called a simple harmonic oscillator, and it undergoes simple harmonic motion: sinusoidal oscillations about the equilibrium point, with a constant amplitude and a constant frequency (which does not depend on the amplitude). If a frictional force (damping) proportional to the velocity is also present, the harmonic oscillator is described as a damped oscillator.  
- *Van Der Pol Oscillator:* In dynamics, the Van der Pol oscillator is a non-conservative oscillator with non-linear damping. It evolves in time according to a second-order ordinary differential equation.   
- *Quantum Harmonic Oscillator:* The quantum harmonic oscillator is the quantum-mechanical analog of the classical harmonic oscillator. Because an arbitrary smooth potential can usually be approximated as a harmonic potential at the vicinity of a stable equilibrium point, it is one of the most important model systems in quantum mechanics. Furthermore, it is one of the few quantum-mechanical systems for which an exact, analytical solution is known.   
To do this, we will use MATLAB to numerically simulate this.    

-------------------------------------------------------------

***Computational Analysis of Coffee Burns***: The modelling of heat related phenomena such as bio-heat transfer is of great importance for development of biological and biomedical technologies. For example, human skin burns, also known as thermal burns, are skin injiries casued by heat, chemicals, electricity or being exposed to radiation. Depending upon the condition and duration of exposure thermal burn may cause sever skin damage. The development of mathematical models has greatly enhanced the ability to analyse various types of bio-heat transfer processes. This surely enable researchers and practitioners in medicine and biomedical sectors to better understand the behaviour of heat transfer affecting the human skin and the proper methods to use to prevent permanent damages. The skin has three main layers: the epidermis, the dermis and the hypodermis. These layers are on top of the muscle and bone. Skin burns can be classified as follows:   
1- First degree (superficial) burns: These burns affect only the epdermis. The burn side appears red and dry with no blisters and is midly painful.   
2- Second degree (partial thickness) burns: These burns involve the epidermis and portions of the dermis. The burn side is red and moist and may be blistered, sollwen and very painul.   
3- Third degree (ful thickness) burns: These burns extend through the dermis and into the hypodermis. The burns side appears patchy and colour ranging from white to brown with a dry lethery texture. Becasue the burn is so deep, it casues little or no pain.  
4- Fourth degree burns: These burns involve the destruction of all layers of the skin. These burns are brown, dry and almost always painless.   
Skin burns could be approximated to a simple 1D heat transfer method if the exposed surface is very large as compared to the burn penetration depth. This allows for studying heat flow in one dimension, i.e., dependent on the depth only. This project aims to solve the 1D heat transfer equation analytically with non-homogeneous boundary conditions and then interpret the behaviour of a bio-heat transfer by sovlving it numerically using Forward Euler Method by writing a code in MATLAB.   

-------------------------------------------------------------

***Chaos in Ordinary Differential Equations***

-------------------------------------------------------------

***PCA and Spring-Mass Systems***: Here 

-------------------------------------------------------------

***Finite Element Analysis of a Flow Over an Airfoil***: The airplane generates lift using its wings. The cross-sectional shape of the wing is called an airfoil.  

-------------------------------------------------------------



