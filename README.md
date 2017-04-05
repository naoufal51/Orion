![ScreenShot](https://cloud.githubusercontent.com/assets/15954923/21268015/28638b76-c3ad-11e6-89d5-6096e56f924a.png)
# Orion
In order to run the notebook showcasing the python version (partial) of the code used for the experiments, please use mybinder or clone the repository and use jupyter on you local machine for launching the notebook present in the code folder. If you choose to use mybinder, you should copy and paste the repository link into the appropriate field present in home page of the website: http://mybinder.org.

We also provide code for CSI measurement using the experiment orchestration tool nepi-ng. It allows configuring the R2lab wireless nodes, run a transmission scenario between a transmitter and a receiver, and collect the CSI logs for angle of arrival (AoA) estimation. This code is tunable for more sophisticated transmission scenarios, for instance the usage of spatial multiplexing. https://github.com/parmentelat/r2lab/tree/public/demos. We also propose an explanatory video for running an AoA experiment scenario https://r2lab.inria.fr/tuto-900-youtube.md#AOA.

# ORION Orientation Estimation Using Commodity Wi-Fi

With MIMO, Wi-Fi led the way to the adoption of antenna array signal processing techniques for fine-grained localization using commodity hardware. MIMO techniques, previously exclusive to specific domains of applications for instance radar systems, allow to consider estimating the orientation in space of a device equipped with COTS Wi-Fi chip. In fact, the availability channel state information measurement in recent Wi-Fi chips, makes it a candidate readily available for experimentation. Accordingly, we propose the ORION system to estimate the orientation (heading and yaw) of a MIMO Wi-Fi equipped object, relying on a joint estimation of the angle of arrival and the angle of departure.

<img src="/code/img/exp_setup2.jpg" alt="Drawing" height="400px">

## Implementation 

### Experimentation Setup
We have implemented the system using the Intel Wi-Fi Link 5300 NICs. Indeed, we have used the 533AN_MMW (Full) model with the following SPS references: 480986-0001 0E, 506679-001 0B, 480986-0010A. The firmware of this COTS Wi- Fi card was modified in order to extract the CSI matrices through the Intel CSI tool for 802.11n HT packets. In case of OFDM systems, as for instance Wi-Fi, we are able to extract a CSI matrix for each sub-carrier. In our case the Wireless NIC offers up to 30 subcarriers. We use the 5 GHz band with 20MHz of bandwidth. We set up our wireless cards in the injection mode, which avoids the need of association with an AP and allows raw Wi-Fi packets transmission.

We setup two uniform linear antenna arrays. We use a classical 5dBi omnidirectional compatible with 2.4 GHz as well as the 5GHz band. As we are using the 5.32 GHz band for our system, we respect a 2.8cm of inter-antenna spacing, which corresponds to half the wavelength. These ULAs are connected to the terminal through extension cables in order to have liberty of movement, which is essential for our experimentation. The cables are 3 meters long and compatible with both the Wi-Fi bands with 2dB signal attenuation for 2.4 GHz band.
 

### Wireless Card Configuration
The transmission is controlled by disabling the antenna selection algorithm and by specifying the desired number of antennas, streams and transmission technique (SM). The rotation is performed using extension cables that connect the coplanar transmitter and receiver antenna arrays to the access points.
In order to specify correctly the antenna mapping as well the number of streams used for the transmission, we change the monitor_tx_rate file, after making sure that all the previous steps for collecting CSI were respected, we could for example launch the following command:

```bash
sudo echo 0x4101 |sudo tee /sys/kernel/debug/ieee80211/phy0/iwlwifi/iwldvm/debug/monitor_tx_rate
```
This command is meant for sending 1 stream on the first antenna. 
In fact the number of streams that are sent is dependent on the MCS index used.

If we want to build a rate index we need to choose and appropriate MCS index corresponding to the number of streams. For example, for sending two streams we choose the MCS index number 9 (make sure that the index is in hexadecimal). For choosing the antenna involved in the transmission we need to apply the appropriate mask:
```
0x04000 for the first antenna
0x08000 for the second antenna
0x10000 for the third antenna
```
We use the HT legacy packet format by using the flag 0x100.

And thus if we want to send 2 streams on the second and third antenna using HT leagacy packets, the appropriate rate would be:
```
0x18109=0x08000+0x10000+0x100+0x9
```
And thus the command used is as follows:
```bash
sudo echo 0x18109 |sudo tee /sys/kernel/debug/ieee80211/phy0/iwlwifi/iwldvm/debug/monitor_tx_rate.

```
For further details check: 
https://github.com/dhalperi/linux-80211n-csitool/blob/csitool-3.13/drivers/net/wireless/iwlwifi/dvm/commands.h#L245-L334

For the transmission we choose the injection mode based on lorcon-old. Lorcon, which stands for Loss Of Radio CONnectivity is library for injecting 802.11 frames https://github.com/dhalperi/lorcon-old.

If you choose to use nepi-ng for running your experiment, you should apply an appropriate rate index in https://github.com/parmentelat/r2lab/blob/public/demos/jobs-angle-measure/angle-measure.sh#L128-L134

## Joint Angle of arrival Angle of Departure Estimation



We jointly estimate AoA and AoD using a 2D MUSIC algorithm on phase corrected CSI matrices collected through out the experiment for a subcarrier.

The figure shown below represent the vaules computed for the pseudo-spectrum. The peak represents the AoA and AoD estimates:

<img src="/code/img/dod-doa.png" alt="Drawing" height="200"/>


```python
import matplotlib.pyplot as plt
%matplotlib inline
from music import *
csi_corr='csi_corr.dat'
csi_target='csi_target.dat'
music(csi_corr, csi_target, 2, 3, 0.5, 0.5,640)
#music_pl(s, 2, 3, 0.5, 0.5,1)
```


![png](/code/img/output_6_0.png)





    array([ -4., -24.])




```python
%matplotlib inline
from music_pl import *
from IPython import display
import time

x = []
y = []
z = []

s = phase_correction(csi_corr, csi_target)    

for i in range(1900):
    x = np.append(x, i)
    o=music_pl(s, 2, 3, 0.5, 0.5,i)
    y = np.append(y, o[0])
    z = np.append(z, o[1])
    plt.gca().cla() 
    plt.plot(x,y,label='DOA')
    plt.plot(x,z,label='DOD')
    plt.ylabel('Angle in degrees')
    plt.xlabel('packet number')
    plt.legend()
    plt.ylim(-90, 90)
    display.clear_output(wait=True)
    display.display(plt.gcf()) 

```


![png](/code/img/output_7_0.png)



![png](/code/img/output_7_1.png)



```python
angle=np.array(np.squeeze([[y],[z]]))
angle.shape
```




    (2, 1942)




```python
import pandas as pd 
df = pd.DataFrame(angle)
df.to_csv("exp_meas3.csv",header=None)
```

###  Smoothing Angle Estimations

As shown in the figure above, the estimates are jittery with a strong presence of outliers. Thus, in order to smoothen our estimates, we apply an outlier detector on the measurement data before using a Kalman filter to mitigate the statistical noise.


```python
plt.figure()
plt.grid(True)
plt.ylabel('Angle in degrees')
plt.xlabel('packet number')
plt.plot(y,'b',linewidth='1',label='DOA') 
plt.plot(z,'r',linewidth='1',label='DOD') 
plt.tick_params(labelsize=14)
plt.legend()
plt.ylim(-90, 90)
plt.show()
```


![png](/code/img/output_12_0.png)



```python
%matplotlib inline
import matplotlib.pyplot as plt
from kalman_fil import *

DOA, DOD,angles1,angles2=kalman_fil("exp_meas.csv")
plt.figure(figsize=(7,7))
plt.grid(True)
plt.ylabel('Angle in degrees')
plt.xlabel('packet number')

plt.subplot(2, 1, 2)
plt.plot(DOA,'b',linewidth='2',label='DOA') 
plt.plot(DOD,'r',linewidth='2',label='DOD')
plt.tick_params(labelsize=14)
plt.legend()
plt.grid(True)
plt.ylabel('Angle in degrees')
plt.xlabel('packet number')
plt.ylim(-60, 60)
plt.subplot(2, 1, 1)
plt.plot(angles1,'b',linewidth='2',label='DOA') 
plt.plot(angles2,'r',linewidth='2',label='DOD')
plt.grid(True)
plt.ylabel('Angle in degrees')
plt.xlabel('packet number')
plt.tick_params(labelsize=14)
plt.legend()
plt.ylim(-60, 60)
plt.show()
```


![png](/code/img/output_13_0.png)


So after applying an outlier removal using a Hampel identifier and a Kalman filter on our estimated data we have significantly decreased the noise and increased the overall accuracy.
