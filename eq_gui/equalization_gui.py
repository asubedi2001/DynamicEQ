import tkinter as tk
import numpy as np
import pyaudio
import wave
import os

p = pyaudio.PyAudio()
# playing -> is audio stream currently being output to user
cur_playing = False
paused = True
audio_file_select = None
audio_stream = None
wf = None

def pause_play():
    global cur_playing, paused, audio_stream, wf
    print(f'cur_playing: {cur_playing}, paused: {paused}, audio_stream: {audio_stream}, wf: {wf}')

    # if audio file is currently being played (paused or unpaused)
    if cur_playing:
        if paused:
            # audio started & was paused, continue audio
            paused = False
            if audio_stream:
                audio_stream.start_stream()
        else:
            # audio started & no paused, pause audio (save where it is)
            paused = True
            if audio_stream:
                audio_stream.stop_stream()
    else:
        # play new audio from beginning
        cur_playing = True
        paused = False

        # load wav file
        filepath = f'./audio/{audio_file_select.get()}'
        if wf:
            wf.close()  # if different wav file opened, close it before assigning new one
        wf = wave.open(filepath, 'rb')  # open new wav file
        
        # open audio stream from wav file
        if not audio_stream:
            audio_stream = p.open(
                format=p.get_format_from_width(wf.getsampwidth()),
                channels=wf.getnchannels(),
                rate=wf.getframerate(),
                output=True,
                stream_callback=stream_callback(wf)
            )
        audio_stream.start_stream()

def stop():
    global cur_playing, paused, audio_stream, wf
    if audio_stream:
        audio_stream.close()
    if wf:
        wf.close()
    audio_stream = None
    wf = None
    cur_playing = False
    paused = True

def stream_callback(wf):
    def callback(in_data, frame_count, time_info, status):
        global cur_playing, paused, audio_stream
        data = wf.readframes(frame_count)
        print(f"Stream callback called. Data length: {len(data)}")
        if wf.tell() >= wf.getnframes():
            print("Playback finished.")
            cur_playing = False
            paused = True
            if audio_stream:
                audio_stream.close()
                audio_stream = None
            wf.close()
            return None, pyaudio.paComplete
        return data, pyaudio.paContinue
    return callback

def equalize():
    pass

# helper for create_gui dropdown, lists songs in ./audio/ folder
def read_audiofiles():
    filepath = './audio/'
    wav_files = [f for f in os.listdir(filepath) if f.endswith(".wav")]
    if wav_files:
        audio_file_select.set(wav_files[0])
        dropdown['menu'].delete(0, 'end')
        for file in wav_files:
            dropdown['menu'].add_command(label=file, command=lambda value=file: audio_file_select.set(value))
    else:
        audio_file_select.set("None")
        dropdown['menu'].delete(0, 'end')
        dropdown['menu'].add_command(label="No files found", command=lambda: None)

# create gui with 3 sliders (low,mid,high), pause/play/stop buttons, and apply eq box
def create_gui():
    window = tk.Tk()
    window.title('Audio Playback GUI for DynamicEQ')
    window.geometry('400x300')

    # create three sliders, representative of bands of frequencies that will be affected
    # by DynamicEQ in the future.
    bass_slider = tk.Scale(window, from_=1, to=10, orient=tk.HORIZONTAL, label="Bass")
    bass_slider.pack()

    mid_slider = tk.Scale(window, from_=1, to=10, orient=tk.HORIZONTAL, label="Mid")
    mid_slider.pack()

    treble_slider = tk.Scale(window, from_=1, to=10, orient=tk.HORIZONTAL, label="Treble")
    treble_slider.pack()

    # create dropdown for .wav file selection
    global audio_file_select, dropdown
    audio_file_select = tk.StringVar()
    dropdown_frame = tk.Frame(window)
    dropdown_frame.pack(pady=10)
    
    dropdown_label = tk.Label(dropdown_frame, text="Select a Song:")
    dropdown_label.pack(side=tk.LEFT)

    dropdown = tk.OptionMenu(dropdown_frame, audio_file_select, [])
    dropdown.pack(side=tk.LEFT)
    read_audiofiles()  # Populate the dropdown on startup

    # create play and pause buttons in button frames
    button_frame = tk.Frame(window)
    button_frame.pack(pady=20)

    play_button = tk.Button(button_frame, text="Play/Pause", command=pause_play)
    play_button.pack(side=tk.LEFT, padx=10)

    # create stop button
    stop_button = tk.Button(button_frame, text="Stop", command=stop)
    stop_button.pack(side=tk.LEFT, padx=10)

    window.mainloop()

if __name__ == '__main__':
    create_gui()