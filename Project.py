__author__ = 'Nathan Davis and Brian Bickford'

from os import listdir
import matplotlib.pyplot as plt
import matplotlib.path as path
from scipy.stats import linregress
from PIL import Image
import math
import copy



def main():
    width, height = 401, 401
    folders = listdir("data")
    veg_counts = []
    urb_counts = []
    wat_counts = []
    for folder in folders:
        # Get position
        # Position for 1999: (3299, 3299)
        # Position for 2009: (3069, 3219)
        x = int(input("What is the X position for " + folder + "?"))
        y = int(input("What is the Y position for " + folder + "?"))
        pos = (x, y)

        # Get reflectance data from images - comment this line if using text files
        bands = find_bands_from_meta("data/" + folder, pos, width, height)



        # Calculate results of Tasseled Cap transformation
        print("Performing Tasseled Cap transformation...")
        brightness, greenness, wetness = tasseled_cap(bands)

        # If the folder name is "1999", compare results to ENVI's
        if folder == "1999" or folder == "2009":
            # Read reflectance data from text files - comment this line if using image files
            text_bands = read_from_text("ScaledReflectance" + folder + ".txt")
            text_tasseled_cap = read_from_text("tasseled" + folder + "data.txt")

            print()
            print("Root mean square deviations for " + folder + ":")
            print("Reflectance Bands:")
            colors = ['b', 'g', 'r', 'c', 'm', 'y']
            for i in range(len(bands)):
                print(root_mean_square(bands[i], text_bands[i]))
                plotlines(bands[i], text_bands[i], colors[i])
            plt.title("Reflectance Comparison: " + folder)
            plt.xlabel("Our Results")
            plt.ylabel("ENVI's Results")
            plt.axis([0, 1, 0, 1])
            plt.savefig("Bandscomparison" + folder + ".png", format='png')
            plt.show()
            print("Tassled Cap Results:")
            print(root_mean_square(brightness, text_tasseled_cap[0]))
            plotlines(brightness, text_tasseled_cap[0], 'r')
            print(root_mean_square(greenness, text_tasseled_cap[1]))
            plotlines(greenness, text_tasseled_cap[1], 'g')
            print(root_mean_square(wetness, text_tasseled_cap[2]))
            plotlines(wetness, text_tasseled_cap[2], 'b')
            print()
            plt.title("Tasseled Cap Comparison: " + folder)
            plt.xlabel("Our Results")
            plt.ylabel("ENVI's Results")
            plt.axis([0, 1, 0, 1])
            plt.savefig("Tasseledcomparison" + folder + ".png", format='png')
            plt.show()

        # Select vegetation, urban, and water pixels from a Tasseled Cap scatter plot
        vegetation = selectpoints(brightness, greenness, "Vegetation")
        urban = selectpoints(brightness, greenness, "Urban")
        water = selectpoints(brightness, greenness, "Water")

        #Assign colors to each point for the final plot and map
        colors = []
        vegetation_count, urban_count, water_count = 0, 0, 0
        for value in list(zip(vegetation, urban, water)):
            if value[0]:
                vegetation_count += 1
                colors.append('g')
            elif value[1]:
                urban_count += 1
                colors.append('r')
            elif value[2]:
                water_count += 1
                colors.append('b')
            else:
                colors.append('y')
        if vegetation_count == 0:
            print("No vegetation cover selected, ending processing.")
            return
        if urban_count == 0:
            print("No urban cover selected, ending processing.")
            return

        # Show a scatter plot with selected points colored according to their type
        print("Graphing selected points...")
        graph2D(brightness, greenness, colors)
        pixels = len(vegetation)

        # Calculate percent of pixels belonging to each category
        vegetation_percent = 100 * (vegetation_count / pixels)
        urban_percent = 100 * (urban_count / pixels)
        water_percent = 100 * (water_count / pixels)
        veg_counts.append(vegetation_percent)
        urb_counts.append(urban_percent)
        wat_counts.append(water_percent)
        print()
        print("Results for " + folder + ":")
        print("Total pixels: " + str(pixels))
        print("Pixels of vegetation: " + str(vegetation_count))
        print("Percent vegetation: %{0:.2f}".format(vegetation_percent))
        print("Pixels of urban: " + str(urban_count))
        print("Percent urban: %{0:.2f}".format(urban_percent))
        print("Pixels of water: " + str(water_count))
        print("Percent water: %{0:.2f}".format(water_percent))
        print()

        # Show map of region with selections highlighted
        show_selected_map(colors, width, height)

    # Calculate changes in vegetation and urban coverage between each year
    years = []
    for folder in folders:
        years.append(int(folder))
    for i in range(len(years)-1):
        veg_change = 100 * veg_counts[i+1] / veg_counts[i]
        print("Vegetation increased by %{:.2f} between {} and {}."
              .format(veg_change - 100 if veg_change >= 100 else 100 - veg_change, years[i], years[i+1]))
        urb_change = 100 * urb_counts[i+1] / urb_counts[i]
        print("Urban cover decreased by %{:.2f} between {} and {}."
              .format(urb_change - 100 if urb_change >= 100 else 100 - urb_change, years[i], years[i+1]))



# Read reflectance values from each of the reflectance text files produced by ENVI
def read_from_text(filename):
    # Open file
    file = open(filename, mode='r', encoding="utf-8", newline='')
    lines = file.readlines()
    if lines == []:
        print("Could not read file.")
        return
    points = []
    bands = []
    # Add points from the file
    for line in lines:
        if line.startswith(";"):
            continue
        line = line.split()
        if line == []:
            bands.append(points)
            points = []
            continue
        for num in line:
            points.append(float(num))
    return bands


def root_mean_square(list_a, list_b):
    pairs = list(zip(list_a, list_b))
    sum = 0
    for pair in pairs:
        sum += (pair[0] - pair[1])**2
    return (sum / len(pairs))**(0.5)



# Use the metadata file to locate and read the Landsat images
def find_bands_from_meta(folder, pos, width, height):
    # Find meta information file
    meta_file = ""
    for file in listdir(folder):
        if file.endswith("MTL.txt"):
            meta_file = folder + "/" + file
    # Read band file information from meta file
    band_info = read_meta_file(meta_file)

    # Get reflectance data from band files
    print("Calculating Reflectance...")
    bands = []
    for band in band_info:
        bands.append(reflect(folder + "/" + band[0], band[1], band[2], band[3], pos, width, height))
    return bands


# Read the raw values from a Landsat image file
def read_meta_file(filename):
    # Open metadata file
    file = open(filename, mode='r', encoding="utf-8", newline='')
    lines = file.readlines()
    band_files = []
    radiance_mult = []
    radiance_add = []
    # Add find image file names
    for line in lines:
        words = line.split()
        if words[0].startswith("FILE_NAME_BAND_"):
            band_files.append(words[2].strip("\""))
        if words[0].startswith("RADIANCE_MULT_BAND_"):
            radiance_mult.append(float(words[2]))
        if words[0].startswith("RADIANCE_MINIMUM_BAND_"):
            radiance_add.append(float(words[2]))
    esun = [1957, 1826, 1554, 1036, 215, 0, 80.67] # from handbook

    # Zip together the information to interpret each band
    band_info = list(zip(band_files, radiance_mult, radiance_add, esun))
    band_info.pop(5)    # Remove band 6
    return band_info


# Perform the Tasseled Cap calculations
def tasseled_cap(bands):
    # Reference: Introductory Digital Image Processing: A Remote Sensing Perspective, Third Edition, John R. Jensen
    brightness = (0.2909, 0.2493, 0.4806, 0.5568, 0.4438, 0.1706)
    greenness = (-0.2728, -0.2174, -0.5508, 0.7221, 0.0733, -0.1648)
    wetness = (0.1446, 0.1761, 0.3322, 0.3396, -0.6210, -0.4186)

    # Create lists with length equal to the size of the bands
    brightness_values = [0 for i in range(len(bands[0]))]
    greenness_values = [0 for i in range(len(bands[0]))]
    third_values = [0 for i in range(len(bands[0]))]

    # Calculate brightness, greenness, and wetness for each pixel
    for j in range(len(bands)):  # Loop through each band
        for i in range(len(bands[0])):  # Loop through each point
            # Multiply point i in Band j by each coefficient and add to total value for that point
            brightness_values[i] += bands[j][i] * brightness[j]
            greenness_values[i] += bands[j][i] * greenness[j]
            third_values[i] += bands[j][i] * wetness[j]
    return (brightness_values, greenness_values, third_values)


# Make a graph with 2 axes
def graph2D(x, y, colors = []):
    if colors != []:
        plt.scatter(x, y, color=colors, marker='.')
    else:
        plt.scatter(x, y, marker='.')
    plt.xlabel("Brightness")
    plt.ylabel("Greenness")
    plt.axis([0.07, 1.1, -0.02, 0.38])
    plt.savefig("2D Plot.png", format='png')
    plt.show()


# Make a graph with 2 axes
def plotlines(x, y, color):
    line = linregress(x, y)
    slope = line[0]
    intercept = line[1]
    plt.plot([0, 1], [intercept, intercept+slope], color)
    plt.scatter(x, y, c=color)
    # plt.axis([0, 1, 0, 1])
    # plt.show()


# Read the raw values from a Landsat image file
def reflect(filename, mult, add, esun, pos, width, height):
    im = Image.open(filename)
    # x1, y1 = (1252, 3240)  # Use these for the text file's region
    x2, y2 = (pos[0] + width, pos[1] + height)
    image = im.crop((pos[0], pos[1], x2, y2))  # get a subset of the image
    flat = list(image.getdata())
    points = []
    for point in flat:
        if point > 0:  # Skip solid black pixels
            l = (mult * point) + add  # Numbers come from the MTL text file under Radiometric rescaling.
            r = (math.pi * l * (1.01429)**2)/(esun*math.cos(math.radians(55.84171719)))
            r = round(r, 4)
            points.append(r)
        else:
            points.append(0)
    return points


def selectpoints(x, y, selection_type):
    plt.scatter(x, y, marker='.')
    plt.xlabel("Brightness")
    plt.ylabel("Greenness")
    plt.title("Select " + selection_type + ":")
    vertices = []
    # Define the event for when a user clicks the figure
    def onclick(event):
        vertices.append([event.xdata, event.ydata])  # add a vertex to the list for the location of the click
        print([event.xdata, event.ydata])  # print the location of the click for the user's reference
    plt.connect("button_press_event", onclick)
    plt.axis([0.07, 1.1, -0.02, 0.38])
    plt.savefig("2D Plot.png", format='png')
    print("Polygon vertices for " + selection_type.lower() + ":")
    plt.show()
    print()
    if vertices == []:
        return []
    # Create a path object from the vertices and select the points within it
    polygon = path.Path(vertices)
    selected = polygon.contains_points(list(zip(x, y)))
    return selected

def show_selected_map(colors, width, height):
    # Build a list with color values for each pixel
    color_values = []
    for color in colors:
        if color == 'r':  # The pixel is urban, use red
            color_values.append((255, 0, 0))
        if color == 'g':  # The pixel is vegetation, use green
            color_values.append((0, 255, 0))
        if color == 'b':  # The pixel is water, use blue
            color_values.append((0, 0, 255))
        if color == 'y':  # The pixel is unselected, use gray
            color_values.append((100, 100, 100))
    # Create a new image with these pixels and display it
    im = Image.new("RGB", (width, height))
    im.putdata(color_values)
    im.show()


if __name__ == "__main__":
    main()
