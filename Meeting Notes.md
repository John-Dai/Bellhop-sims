##Mar 22nd

Vaishakhi is working on the decoder & mapper 
	Amplitude degradation as compared to the original signal 
	Hitting a hard surface – changes in phase and amplitude 
	Pulses locations increase and decrease
	Location of transmitter, but not angle, changes 
	Taking the changes in beams and stitching them together to create a map 	Adding strips of land over the area that we see
	Topographical maps
	Only one transmitter one receiver, close together
She suggests: 
Read the documentation 
Look at the code 
Ask the people who made it for help 
Need someone who has worked on mixed file wrappers to help us? – for a week or two weeks – to help us  - V will talk to martin about this 
3d impulse code is different – wrapper 
2d – purely matlab code 

Audio file is most important 
	Main challenge is reading the data 
Harmonics? 
	Time? 
	V agrees with M’s idea of summing up time delayed 
A variety of 2d 3d bty files 
	2D – easy to make using the drawbty code, just need to make more 
	3D – need to write the code to draw this 
Documentation for how to run bellhop 
What v wants: 
	Sound profile as the pod goes down 
		50 transmitters and 50 receivers – map pairs of receivers & transmitters as they go down – vertical variation 
		Horizontal movement difficult 

V looking at adding noise 

V’s algorithm uses one sender/receiver pair
She wants just a wav or data audio file - considering that .wav is what the prof suggested in the meeting we’ll probably be good 
She’ll give us a link to the algorithm
We should ask hardware team for a sample output


##Mar 21st

Audio file conversion 
	Add noise in the ocean – standard ocean noise of 80db – Professor has a way to convert into volts that he found working with a sound engineer trying to test the range of the sensors for the hardware team 1600m range
	Send out a simple signal 
		Arbitrary magnitude (1) sine wave at 30kHz, for 30 cycles, 
		





Convolve this data with a matching filter – real time 



	

Multiply original signal by relative amplitude and add to original signal (if you’re using sound pressure – you can add it, if we’re in dB we cannot add)
If we multiply by the numbers Bellhop gives us, we need to find the units – what does Bellhop want out of this ‘
	dB references are difficult – divided by what to get dB?
	Kilo Pascals unit for pressure 
May use a chirp for better results 
.wav files – incorporate all the information self contained within the wav file, more transportable 

Paper:
	Prof. vastly prefers Google Docs to Latex – not a source document, unless we have a compelling reason  
	A chance for the prof. to give us feedback for the final paper 
	Talking about the paper, making plans
	Main goal: helping Vaishakhi 
	Important references: reference any tools you use, like MATLAB or bellhop, reference everything we use to think about the problems, lots of references to valid sources is important, citing links in Google Drive pages to make sure links don’t break 

	Pick cool images to put in the document 
	
Competing this fall 



##Feb 28th

John Dai made a beginning pass at pond bathymetry, modeling an uneven seabed and a wall.

Professor Brooke says we need 360 emanation from our sensors and more rays (see John’s

screenshots on spreadsheet). The goal is to be able to compare our simulations to what the

hardware team gets from the pond.

Ask E if you have questions about the parameters in the bottom half of the simple.env file (found

in Google Drive folder \Bellhop\Spring 2017)

Begging for money (see email):

We’re asking for funding from certain companies that might be interested. This involves

contacting the companies and pitching our XPRIZE entry to them, either in person or by sending

them a video.

Things Professor Brooke needs volunteers for:

Editing the document he sent us

Contacting the companies

Creating slides for a presentation

Presenting our XPRISE entry (irl or in a video)

Recording/editing the video if we make one

Audio file progress:

##Feb 21st

John Dai got a reflecting bottom 
need to be changing parameters for ENV files and understanding impulse responses 
V's code right now finds out locations based on time series data 
V wants us to make time series data 
multiple receivers at different rangers and depths - only want one receiver closest to the transmitter 
combine time series data: 
sound files 
need to understand output 
converting .arr into audio file (data in time) 
monotone sine waves at 28/30kHz
Can do any receiver 
Is bellhop supposed to support sound files? 
sound file extension should be 
not big 
not compressed  
Can we control how long the sound is? Short chirp or long pulse or continuous sound?  
need to generate large simulations to test V's code 
testing moving sensors and receivers - multiple ways to do it which V is testing 
cool idea: sending from one location and receiving at different depths and ranges to simulate movement 
2D (moving to 3D later) 
bottom features characterization
organize in a spreadsheet by features & parameters 
don't change: transmission frequency
note what you change and how 
splitting teams: 
2D .arr to audio file: Henry, Chris
bottom features characterization: Eeyi, John Dai 

For later, but important: 
2D to 3D 
movement 


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


