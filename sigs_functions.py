import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy.signal import freqz
from scipy.io import wavfile
from scipy.io.wavfile import write, read
import wave
from dataclasses import dataclass

@dataclass
class SignalInfoOutput:
    data: list
    samplerate: int
    duration: float|int





def signal_generator(
        signal_type: str | None, #programming mate thought me to define my variable types, making it easier when implemnting the functions in other document 
        sampling_rate: int | float,
        duration, 
        frequency,
        amplitude,
        phase
    ) -> np.ndarray:
    signal_type = signal_type if not None else 'sinusoid'
    if frequency < 0:
        raise ValueError("frequency below 0")
    

    t = np.linspace(0, duration, int(sampling_rate * duration), endpoint=False)
    #create sinusoid
    sinusoid = amplitude * np.cos(frequency * 2*np.pi * t + phase)

    #create square wave 
    square_wave = []
    for number in sinusoid:
        if number > 0:
            square_wave.append(amplitude)

        elif number < 0:
            square_wave.append(-amplitude)

        elif number == 0:
            square_wave.append(0)

    #create triangel
    triangel_wave = []
    dt = t[1] - t[0] #change in time, is constant, so only one dt needed
    integral_square = 0

    for value in square_wave:
        integral_square += value*dt #value*dt is the integral of the square wave, which is a triangle wave
        triangel_wave.append(integral_square) 
        
    #create saw tooth
    sawtooth_wave = sp.signal.sawtooth(2*np.pi * frequency * t)


    if signal_type == 'sinusoid':
        return sinusoid
    elif signal_type == 'square wave':
        return square_wave
    elif signal_type == 'triangel wave':
        return np.array(triangel_wave)
    elif signal_type == 'saw tooth':
        return np.array(sawtooth_wave)
    else:
        raise ValueError("Wrong type hilsen sigs")




def signal__to_plot(signal, sampling_rate, duration, frequency, amplitude, phase, x_lim, signal_type):
    t = np.linspace(0, duration, int(sampling_rate * duration), endpoint=False)
    plt.figure()
    plt.xlim(0, x_lim)
    plt.plot(t, signal, label = signal_type)
    #plt.stem(t, signal_type)
    plt.legend()
    plt.show()



#Co-pilot suggestsed using masks to create a piecewise function. Boolean masks will allow us to look at the t data set and change it as we work though our t values.
def ADSR_envelope_func(sampling_rate,
                        duration, attack_t,
                        decay_t, sustain_t,
                        release_t, 
                        sustain_level): #GG! parameters x_t should add up to duration
    
    t = np.linspace(0, duration, int(sampling_rate * duration), endpoint=False)


    ADSR_envelope = np.zeros_like(t) #make a new enevlope array with same length as t, full of zeroes. 

    #defininfg attack section
    attack_mask = t < attack_t #so attack this prositive growing line goes from 0 to 0.1, so we have defined how long our attack is
    #tilt are based on linear ramps, t/duration, just like in a linear function. 
    attack_tilt = t[attack_mask]/attack_t
    ADSR_envelope[attack_mask] = 1 * attack_tilt #all values that are true will then get a new value of t[attack_maks]/attack time, this is t-vlaue/0.1,  we have normalized our time_axis, so we end with an amplitude of 1. 

    #definig decay section
    decay_mask = ((t >= attack_t) & (t < attack_t + decay_t)) #setiting the limits fot the decay time
    decay_tilt = ((t[decay_mask] - attack_t)/ decay_t ) #subtarct with attack_t to make it start at zero
    ADSR_envelope[decay_mask] = (1 - (1 - sustain_level) * decay_tilt) # we want decay to start at 0 relative to the start of decay section, so we subtrack t[decay_mask- attack time. dividing by decay time, we normalize decay section. We say 1 - sustain level, as this is how far we need to drop from the start of the decay to the start of sustain level. Lastly we subtract it all from one, so that the array starts at where the attck section ended. 

    #Defining sustain section
    sustain_mask = ((t >= attack_t + decay_t) & (t <  attack_t + decay_t + sustain_t))
    ADSR_envelope[sustain_mask] =  sustain_level 

    #defining release sectionr
    release_mask = ((t >= attack_t + decay_t + sustain_t ) & (t < attack_t + decay_t + sustain_t + release_t))
    release_tilt = ((t[release_mask] - (attack_t + decay_t + sustain_t))/release_t)

    ADSR_envelope[release_mask] = sustain_level * (1 - release_tilt) #release is decaying

    return ADSR_envelope


def filter_element(filter_type, order, sampling_rate, frequency, signal, band_width, cut_off_n):
    """
    returns: signal, a, b 
    """


    cut_off = lambda n: frequency*2**(n/12)
   # cut_off_norm = cut_off/(sampling_rate/2) #nomralized to nyquist, sampling_rate/2 = nyquist, the sp.butter want it normalized compared to the nyquist freq. (apparently you don't need to normalize to nyquist)

    print(cut_off(cut_off_n))
    cut_off_2 = cut_off(cut_off_n) + band_width #nomralized to nyquist

    


#butterworth filter design
    lb, la = sp.signal.butter(order, cut_off(cut_off_n), btype = 'lowpass', fs = sampling_rate) #when i give sp.butter fs =, it will automatically normalize my cut_off to nyquist frequency
    hb, ha = sp.signal.butter(order, cut_off(cut_off_n), btype = 'highpass', fs = sampling_rate)
    bp_b, bp_a = sp.signal.butter(order, [cut_off(cut_off_n), cut_off_2], btype = 'bandpass', fs = sampling_rate)

    #pass signal through  our filter, sosfilter is generally more stable than lfilter for higher order
    high_signal = sp.signal.lfilter(hb, ha, signal)
    low_signal = sp.signal.lfilter(lb, la, signal)
    band_signal = sp.signal.lfilter(bp_b, bp_a, signal)
    
    #return wanted signal
    if filter_type == 'lowpass':
        return np.array(low_signal), la, lb
    
    elif filter_type == 'highpass':
        return np.array(high_signal), ha, hb
    
    elif filter_type == 'bandpass':
        return np.array(band_signal), bp_a, bp_b
    


#From hands on 6, where we plotted the frequency, resuing parts of code
def plot_frequency_response(a,
                            b,
                             sampling_rate, 
                             x_lim):
    freq_response = sp.signal.freqz(b, a, fs = sampling_rate)
    wf, wh = freq_response
    freq_response = sp.signal.freqz(b, a, fs=sampling_rate)
    wf, wh = freq_response

    #print(np.max(np.abs(wh))) #uncomment to check if it is working, rough debug
    #print(np.min(np.abs(wh)))
    #print(f'a= {a}')
    #print(f'b= {b}')
    wh_dB =  20*np.log10(np.abs(wh) + 1e-10 )#we want to plot in decibel, this makes more sense for our purposes as we want to make sounds. 
    plt.figure()
    plt.xlim(0, x_lim) #using an x-lim, to get a better view of the filter during the cut off, however the x_lim will depend on esch function you plot. 
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("dB")
    plt.plot(wf, wh_dB, label = 'Frequency Response') 
    plt.legend()
    plt.show()


def allpass_filter(signal, alpha: float):
    """
    returns: filtered_signal, a, b
    """
    b = np.array([alpha, 1]) #defining numerator as alpha + x^-1
    a = np.array([1, alpha]) #defining denomenator as 1 + aplha*z^-1
    alpass_signal = sp.signal.lfilter(b, a, signal) #using lfilter to pass signal through all pass, using the correct numerator and denomerator. 

    return alpass_signal, a, b


#plotting function for phases asked co-pilot on how to plot in phase spectrum. 
def plot_phase_response(a, b, sampling_rate):

    w, h = sp.signal.freqz(b, a, fs=sampling_rate)

    phase = np.unwrap(np.angle(h))

    plt.figure()
    plt.plot(w, phase)
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Phase (radians)")
    plt.title("Phase Response")
    plt.grid()
    plt.show()


def reverb(signal, sampling_rate):
    #choosing to add 4 delays around 20 ms apart, finding the integer at what point to start the delay in the original signal
    delay_1 = int(0.02 * sampling_rate) #first 20ms delay
    delay_2 = int(0.04 * sampling_rate) # 40ms delay
    delay_3 = int(0.06 * sampling_rate) # 60ms delay
    delay_4 = int(0.08 * sampling_rate) # 80ms delay

    reverbed_signal = signal.copy()

    reverbed_signal[delay_1:] += 0.8 * signal[0 : len(signal) - delay_1] # we add a delay to the signal, and lower the amplitude, so it sounds more like an echo
    reverbed_signal[delay_2:] += 0.6 * signal[0 : len(signal) - delay_2] 
    reverbed_signal[delay_3:] += 0.4 * signal[0 : len(signal) - delay_3] 
    reverbed_signal[delay_4:] += 0.2 * signal[0 : len(signal) - delay_4] 

    #Now we have a "clear" echo of the signal but to get a better effect we can use the All-pass filter to "smudge" the sounds together.
    #we pass the revebed signal several times through the all pass filter changing the alpha
    #I decided to just keep the alpha  as set values as I think adding these as parameters ->
    #-> in the function was complicating things unnecessarily, and it is simply faster and simpler to use the reverb function as this. 
    reverbed_signal = allpass_filter(reverbed_signal, 0.7)[0]
    reverbed_signal  = allpass_filter(reverbed_signal, 0.6)[0]
    reverbed_signal  = allpass_filter(reverbed_signal, 0.5)[0]

    print(type(reverbed_signal))

    return reverbed_signal


#co-pilot save_signal function, to make the synthesizing process simpler
def save_signal(signal, sampling_rate, filename: str):
    signal = signal / np.max(np.abs(signal))  # signal values are normalized to [-1,1] so they fit into the  [-32768 ... 32767] wav file's range/ 
    audio = np.int16(signal * 32767) #here we see the normalized values get resized for the wav file, this will also prevent clipping when we convert the files. 

    write(filename, sampling_rate, audio)




#from hands_on_1 exercise 2, I used a similar function to import my data, and get relevant paramenters
def signal_info_data(file_path): 
    samplerate, data = wavfile.read(file_path)
    number_of_channels = 1 if data.ndim == 1 else data.shape[1]

    print(f'number of channels: {number_of_channels}') #here we check if the data is from stereo, if it a has (x, 2) it means there are two columns, 2 audio channels, one for each ear, then we can handle this later
    print(f'samplerate = {samplerate}, and number of samples = {len(data)}')#so song has sample rate of 44100, so it changes frequency 44100 per second

    #duration is equal to the length of the sound file, we take all the samples and divide by sample per second, to get the amount of seconds
    duration = len(data)/samplerate 

    print(f'the duration of the sound file is: {duration} s')

    #use the open function in wave library which gives us the numbet of bytes, times by 8 to get the bit size.
    bit_size_per_sample = (wave.open(file_path, mode=None).getsampwidth())*8

    print(f'bit size per sample: {bit_size_per_sample}') 
    out = SignalInfoOutput(data = data, samplerate=samplerate, duration=duration)
    return out

