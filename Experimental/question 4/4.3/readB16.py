
import os
import struct
import numpy as np

def readB16(dirPath, filename=None):

    # Parse optional input
    if filename is not None:
        filename = os.path.join(dirPath, filename)
    else:
        filename = dirPath

    # Open the file
    try:
        f = open(filename, "rb")
    except OSError:
        raise FileNotFoundError(f"Could not open file: {filename}")

    with f:
        # Check that it is a PCO file
        filetype = f.read(4)
        if filetype != b"PCO-":
            try:
                ft = filetype.decode("ascii", errors="replace")
            except Exception:
                ft = str(filetype)
            raise ValueError(f"Wrong filetype: {ft}")

        # Get image dimensions:
        fileSize, headLength, imgWidth, imgHeight = struct.unpack("<4i", f.read(16))

        # Look into the extended header, throw error if color image
        extHeader = struct.unpack("<i", f.read(4))[0]
        if extHeader == -1:
            colorMode = struct.unpack("<i", f.read(4))[0]
            if colorMode != 0:
                raise ValueError(
                    "Color image detected. Only b/w images have been tested with this function"
                )

        # Get the image (match MATLAB column-major fread + transpose)
        f.seek(headLength, os.SEEK_SET)

        # Read exactly width*height pixels as uint16
        n = imgWidth * imgHeight
        data = np.fromfile(f, dtype=np.dtype("<u2"), count=n)
        if data.size != n:
            raise EOFError("File ended before all pixels were read.")

        # MATLAB: fread(...,[imgWidth,imgHeight],'uint16')'  -> (imgHeight x imgWidth)
        image = data.reshape((imgWidth, imgHeight), order="F").T

    return image