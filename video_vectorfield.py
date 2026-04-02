import matplotlib.pyplot as plt
import random
import numpy as np
import os
import cv2
import glob

def get_frame_size(image_path):
    """ Reads an image and calculates
    the width and length of the images,
    which will be returned """
    frame = cv2.imread(image_path)
    height, width, layers = frame.shape
    frame_size = (width, height)
    return frame_size


def video_from_images(folder='.', 
                      video_name = 'video.avi',
                      suffix='png',
                      prefix='pic_',
                      reverse=False,
                      length_of_video_in_seconds=None,
                      codec = cv2.VideoWriter_fourcc(*'DIVX')):
    """ The function creates a video from all the images with
    the suffix (default is 'png' and the prefix (default is 'pic_'
    in the folder 'folder'. If 'length_of_video_in_seconds' is set
    to None, it will be the number of images in seconds. If a positive
    value is given this will be the length of the video in seconds.
    The function assumes that the the shape of the first image is
    the one for all the images. If not a warning will be printed
    and the size will be adapted accordingly! """
    
    images = []
    pattern = f'{folder}/{prefix}*{suffix}'
    for fname in glob.glob(pattern):    
        images.append(fname)
    if (len(images) == 0):
        print(f"No images found with pattern {pattern}")
        return

    images.sort(reverse=reverse)
    if length_of_video_in_seconds is None:
        # each image will be shown for one second
        length_of_video_in_seconds = len(images)
    
    # calculate number of frames per seconds:
    frames_per_second = len(images) / length_of_video_in_seconds
    frame_size = get_frame_size(images[0])

    video = cv2.VideoWriter(video_name, 
                            codec, 
                            frames_per_second, 
                            frame_size)  
    
    for image in images:
        im = cv2.imread(image)
        height, width, layers = im.shape
        if (width, height) != frame_size:
            print(f'Warning: {image} resized from {(width, height)} to {frame_size}')
            im = cv2.resize(im, frame_size)
        video.write(im)
    cv2.destroyAllWindows()
    video.release()

video_from_images(prefix='plot', suffix="jpg", length_of_video_in_seconds=10)
