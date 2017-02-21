##Feb 14th

Having gotten Vaishakhi’s code working, we’re interpreting some of the results from the .arr graph. 
3 graphs:
•	fixed range
•	fixed depth across a range (like looking at attenuation)
o	continuous? 
Initially doing 30 beams, running 300 beams – 
	How do they decide where the beams go? 
	John Dai’s graph models a beam stuck out in free space – not free space, but an empty ocean 
	Directivity information: alpha1, 2; beta 1, 2 a 2d system aimed at a 30 degree angle 
Two modes – two times hitting the sensor at a particular range and time – from surface echo? 
Each dot is a ray hitting the director, 
Not reflected 
Only happens after a while in time 

Understanding the setup for the simulation and then try to create something
Simpler test scenario’
Deflection layer in the top of the ocean  - each depth has a certain speed 
Prof wants a simpler model – with no ocean – making the ocean super deep 
Simpler model: 
-	No doubling (so that must have been from the echoes of the surface and ocean floor)
-	Dying out faster 
-	Dots along a vertical line mark where the receiver got on the Gaussian beam as well as attenuation 
-	Changing across the range is different receivers, and what they get over time 
Add one object to reflect off of 
	Changing the beam to point either on the top or on the bottom 
	To look at the impulse responses – one receiver examine the time evolution 
How many rays does it take to create an interesting distribution 
Slowly increasing the complexity of our models to compare with the hardware team 
Sources and receivers stepping down – to understand movement 
Alpha changes the direction of the collection of rays: 75 90 to shine down 
Plotarr measures the arrivals for one specific receiver 
Looking at time response we can see the flat surface 
1.	Notes
2.	Adding reflection to the bottom
3.	Creating the impulse response from the .arr file 
4.	Adding features to the bottom of the ocean 
		




##Feb 7th
V – converting data into SONAR data, not yet reconstructed into 
A set of simulations and come up with ways they could fail – whale false flags 
Data+ quick analysis 
vaishakhi.mayya@duke.edu  Gross Hall 330 
Our job to create simulations of the environment & creating sonar data 
First step: 
Realistic echo – one pulse out, one pulse back in, comes from the impulse response 
2 months: simulate something, V tries to reconstruct it, sees what match up 
	Time series signal – modulated sound wave 
	For now, ARR files will be good, if we can generate some kind of echo 

How big does the bellhop think the receiver is? Gaussian beams vs rays 
Big sphere as obstacle – more likely to produce ray – curved surface 
If a ray misses when it’s reflected back, do you get a signal? 
Simulating wave propagation with ray tracing and having Gaussian beams that follow it 
Adding noise? 
How many rays are necessary to really sense any details?  - GPU 
Bellhop can’t run on AWS? 
How would you simulate movement? 
How would you simulate the fact that your pod is moving horizontally as well as vertically? 
-	Require change for .env 
Impulse response 2D – leads us to inverse FT find the attenuation, to tell where the stuff is
Send pulse out & get something back 
Code for spotlight radar, strip map radar 

Arr 3D is done by doing 2D cross sections – range in the right direction. 

Team that collected data from the pond & have calibrated sensors


