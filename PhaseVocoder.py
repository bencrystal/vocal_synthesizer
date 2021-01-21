#https://python-sounddevice.readthedocs.io/en/0.3.12/

#Double check Matlab element 1 = Python element 0

#from aetypes import end

import sounddevice as sd
import numpy as np
import math as math
import scipy.io.wavfile as wav

fs = 44100
duration = .1        # seconds
print ("Recording Audio")
myrecording = sd.rec(duration * fs, samplerate = fs, channels = 1,dtype = 'float64')

sd.wait()


print ("Audio recording complete, Play Audio")
sd.play(myrecording, fs)
sd.wait()
print("Play Audio Complete")

#Do I need to include bits, channels, etc?


x = myrecording

n = 1024
hop = n/4
w = n
f = n
# Pitch Modulation

if (np.size(x) > 1):
    x = np.transpose(x)
#end

s = np.size(x)
#s = len(x)


if (np.size(w) == 1):

    if (w == 0):
        w = w + 1
    #end
    halflen = (w-1)/2
    halff = f/2




    halfwin = .5 * (1 + np.cos(math.pi * np.arange(0, halflen)/halflen))               #(0:halflen)/halflen))   AND PI

    #saving all cosin values into array halfwin
    win = np.zeros((1, f),dtype=np.int8)
    acthalflen = min(halff, halflen)


    win = np.ndarray.flatten(win)   #flatten win so elements can be referenced


    #WIN HAS 1024 PLACES BUT CAN'T BE CALLED without being flattened

    #CONVERT HALFF AND ACTHALFLEN TO INTS TO BE REFERENCED BELOW, MAYBE DEFINE EARLIER TO REDUCE NUMBER OF TYPE CHANGES
    #print(acthalflen)
    #[np.arange(0, acthalflen)]


    ######win[np.arange((halff + 1), (halff + acthalflen))] = halfwin[np.arange(0, acthalflen)]

    #win[np.arange((halff + 1), (halff + acthalflen))] = halfwin(np.arange(1, acthalflen))
    ######win[np.arange((halff + 1), (halff - acthalflen+2), -1)] = halfwin[np.arange(0, acthalflen)]




    win[np.arange((int(halff)), (int(halff) + int(acthalflen)))] = halfwin[np.arange(0, int(acthalflen))]

    # win[np.arange((halff + 1), (halff + acthalflen))] = halfwin(np.arange(1, acthalflen))
    win[np.arange((int(halff) + 2), (int(halff) - int(acthalflen) + 2), -1)] = halfwin[np.arange(0, int(acthalflen))]






#end
else:
    win = w
#end

w = len(win);
# now can set default hop
if (hop == 0):
    hop = math.floor(w/2)
#end

c = 1;

# pre-allocate output array
a = 1+np.fix((s-f)/hop)

if a < 1:
    a = 1
#end
#print(a)
#d = np.zeros((1+f/2), 1+np.fix((s-f)/hop))
print(f)
print(type(f))
d = np.zeros((1+f/2, int(a)))
#print(d)

#print(np.size(d))
#print(1+np.fix((s-f)/hop))  #If duration is too small, D won't have enough dimensions (in matlab, between .05 and .1)
#print(s)
#print(f)
#print(hop)

x = np.array(x)
win = np.array(win)


#DO I NEED TO PREALLOCATE U?

for b in range(0, s-f, hop):                              #(b = np.arange(0,s-f,hop)):      #0:hop:(s-f)
    u = win*x(np.arange((b+1), (b+f)))
    t = np.fft.fft(u) #window by window FFTs


    d[:, c] = np.transpose(t(arange(0, 1+f/2)))             #1:(1+f/2))'     #fills in output array
    c = c+1
#end

  # Otherwise, no plot, but return STFT
D = 1.0*d


#print(type(d))
#print((d))
#print(np.size(D))
#print(type(D))
#print(D)

################PVOC
#Calculate the new timebase samples


#print(d)
#print(D)
rows, cols = np.shape(D);



r = .8 #timescale signal r times faster, CHANGEABLE (.8 for major thirds, .6667 for fifths, .5 for octaves
#maybe R3,r3,r5,r8, etc as variable names?

t = np.arange(0,cols-2,r)                               #t = 0:r:(cols-2);






##################NEW PVSAMPLE SELECTION
rows, cols = np.shape(d) #size of stft array
N = 2*(rows-1); #fft size

if (hop == 0):
  # default value
  hop = N/2
#end

# Expected phase advance in each bin
dphi = np.zeros((1, N/2+1))
#print(dphi)
#print((2*math.pi*hop)/(N/np.arange(1,N/2)))
dphi[0,np.arange(2,1+N/2)] = (2*math.pi*hop)/(N/np.arange(1,N/2))



#dphi(2:(1 + N/2)) = (2*pi*hop)./(N./(1:(N/2)));

# Phase accumulator
# Preset to phase of first frame for perfect reconstruction
# in case of 1:1 time scaling

#print(np.angle(d[:,0]))
ph = np.angle(d[:,0])

# Append a 'safety' column on to the end of b to avoid problems
# taking *exactly* the last frame (i.e. 1*b(:,cols)+0*b(:,cols+1))
d = [d,np.zeros((rows,1))]

ocol = 1





#################### PVSAMPLE TIMESHIFTING (for loop trouble)

c = np.zeros((rows, len(t)))
print('c is ')
print(c)
print(t)
tt = 0
i = 0
for i in range(len(t)):    #for (i = np.arange(0,cols-2,r)):

    # Grab the two columns of b

    dcols = d[:, floor(tt) + [1, 2]]
    # dcols = d(:, [1 2]); % Fixes indexing error?

    tf = tt - floor(tt);
    dmag = (1 - tf) * abs(dcols[:, 1]) + tf * (abs(dcols[:, 2]))         #THERE MAY BE ERRORS WITH THE LAST DCOLS
    #calculate phase advance
    dp = angle(dcols[:, 2]) - angle(dcols[:, 1]) - np.transpose(dphi)
    # Reduce to - pi:pi range
    dp = dp - 2 * pi * round(dp / (2 * pi))
    # Save the column
    c[:, ocol] = np.matmul(dmag,math.exp(j*ph))                 #dmag. * exp(j * ph);

    # Cumulate phase, ready for next frame
    ph = ph + np.transpose(dphi) + dp
    ocol = ocol + 1
    tt = t[i]

#end



################# ISTFT
ftsize = n  #second input of istft function
h = hop; #set default hop size for istft

s = np.size(c)
print(s)
if (s[0] != (ftsize/2)+1):
  error('number of rows should be fftsize/2+1')
#end
cols = s(2);





if (length(w) == 1):
    if (w == 0):
    # special case: rectangular window
        win = ones(1,ftsize)
    elif (rem[w, 2] == 0):   # force window to be odd-len
        w = w + 1
    #end
    halflen = (w-1)/2
    halff = ftsize/2
    halfwin = 0.5 * ( 1 + cos( pi * (np.arange(0,halflen)/halflen)))  #Maybe issues with the np.arange, might be single value?
                               #(0:halflen)/halflen));
    win = zeros(1, ftsize);
    acthalflen = min(halff, halflen);
    win[np.arange((halff+1),(halff+acthalflen))] = halfwin(np.arange(1,acthalflen))
    win[np.arange(halff+1,halff-acthalflen+2,-1)] = halfwin(np.arange(1,acthalflen))
    #Make stft-istft loop be identity for 25% hop
    win = 2/3*win;
    end
else:
    win = w;
#end
w = len(win)

# now set default hop
if (h == 0):
  h = floor(w/2);
#end

xlen = ftsize + (cols-1)*h;
x = zeros(1,xlen);

for b in range(0,h*(cols-1),h):

    ft = np.transpose(c[:,1+b/h])
    ft = [ft, conj(ft(np.arange(ftsize/2,2,-1)))]                   #ft = [ft, conj(ft([((ftsize/2)):-1:2]))];
    px = real(ifft(ft));
    x[np.arange(b+1,b+ftsize)] = x(np.arange(b+1,b+ftsize))+np.matmul(px,win)
#end