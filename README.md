**This is an old college project. It was setup for a much older version of OSX and SDL 1.2 and as such, doesen't compile out of the box, but I want to save it on github for prosterity. (I didn't write the rest of this README, my prof did.)**

## Overview
This source code provides a complete example application for viewing amc/asf
motion capture data.

## Controls
The application reads data from a directory tree, expecting one asf file and
possibly several amc files per directory. Once the motion has been loaded,

* PAGEUP/PAGEDOWN selects motion for playback. On the Mac, hold down the fn key and type the up/down arrow keys.
* SPACE to toggle 0x (pause), 0.1x (slow), 1x (normal) speed
* TAB to toggle camera tracking, 
* ESC to quit
* mouse buttons for camera control. One the Mac, the three mouse
buttons are obtained using the normal button, the normal button while
holding down the option key, and the normal button while holding down
the "Apple" key.
* 'd' will dump all the currently loaded motions to associated .global
files. In such a file, each line contains a reference position of the character
on the floor (position.x, position.z) and a reference yaw (position.yaw),
followed by various joint positions in this reference coordinate frame.
(There is a line for each frame of the motion.)

Using the function defined in Vector/Misc.hpp, that means that, for example,
the root position is given:

rotate_by_yaw(root_from_file, position.yaw) + make_vector(position.x, 0.0, position.z)

In other words, to go from the root.xyz in the file (root_from_file) to the
global root position, you'd use the following formulas:

root.x = root_from_file.x * cos(position.yaw) + root_from_file.z * sin(position.yaw) + position.x
root.y = root_from_file.y
root.z = root_from_file.z * cos(position.yaw) + root_from_file.x * -sin(position.yaw) + position.z

## Compiling
This code should compile without problem under linux; the Jamfile setup
under windows is a little trickier, but it should be doable (I use this basic
jamfile for other windows programs). The directory navigation will need to be
rewritten for windows, however. Under OS X, you can use jam, sdl, and libpng
from macports. You can also use Xcode and the OS X version of SDL. 
See compilation directions below.

To compile, use jam (http://www.perforce.com/jam/jam.html). You'll also need
SDL (http://www.libsdl.org/).

## Detailed Notes on Compiling for OS X using Xcode
You will need SDL from http://www.libsdl.org/. In the left column of that
page, find the Download section and click on the SDL 1.2 link. Click on 
the SDL-1.2.13.dmg Runtime Library. A finder window will open with a
virtual disk named SDL. Drag the SDL.framework to /Library/Frameworks.

You should now be able to open amc_viewer.xcodeproj in Xcode to build
and run. Since you may have a newer version of Xcode, the following details
may help you build your own amc_viewer.xcodeproj from "scratch". The 
following should not be necessary if you can use the provided
 amc_viewer.xcodeproj.

You should now be able to run from within Xcode. You can copy the
binary amc_viewer from amc_viewer.app and invoke it from the
command line. These instructions will not produce a 

## Creating an Xcode project from scratch

Next, you should click on the SDL-devel-1.2.13-extras.dmg link to get
templates for Xcode. Drag the templates from the 
SDL-devel-extras/TemplatesForXcode folder to 
/Library/Application Support/Apples/Developer Tools/Project Templates/Application.

Go to Xcode and open a new project. The template to use is 
SDL OpenGL Application. Give the project a name (amc_viewer)
and location, e.g. amc_viewer.r1523. It is 
probably a good idea at this point to actually build and run the default
application before going further. This will assure you that SDL is installed 
and working.

Now, you can modify the project to make the amc_viewer.

Remove the .c and .h files, but do not remove SDLMain.m or frameworks.
Add Browser/*.cpp, Browser/*.hpp, Character/*.cpp, Character/*.hpp, 
Graphics/*.cpp, Graphics/*.hpp, Library/*.cpp, Library/*.hpp,
Vector/*.cpp, and Vector/*.hpp

Now *remove* Texture.cpp and Texture.hpp from the project.
Also, *remove* ReadSkeletonV.cpp, WriteAsfAmc.cpp, WriteAsfAmc.hpp, 
                        WriteBvh.cpp and WriteBvh.hpp from the project.

In Project/Edit Project Settings, add the amc_viewer.r1523 directory
to the Header Search Paths. You want this directory to appear in 
both release and debug configurations, but the default is to edit
just the release version, so be careful.

VectorGL.hpp has type conversion problems, so comment out
the definitions of: 
inline void glVertex(Vector2i const &v, int const &z) and
inline void glVertex(Vector3i const &v)

Now Build and Run the application. Probably it will complain that there
are no motions to browse in the data directory, and the program needs
a font file, so you're not done yet.

Create a data directory. Copy some motions to it, e.g. copy both 
93_06.amc and 93.asf to the directory named data.

In Xcode, go to Project/Edit Active Executable and in the General tap,
Set the working directory to: Custom directory: dist

Then under the Arguments tab, add the argument ../data
(This assumes your directory structure is something like:
      amc_viewer/dist -- run directory, contains gentium.txf
      amc_viewer/data -- contains motion files

Again, be careful to configure both deployment and debug versions
of your program.

You can now run from within Xcode or copy the binary amc_viewer from
amc_viewer.app to dist and run from a terminal window.

## Data
A good source for motion capture data is http://mocap.cs.cmu.edu/.

## Copyright
All source is copyright Jim McCann unless otherwise noted. Feel free to use
in your own projects (please retain some pointer or reference to the original
source in a readme or credits section, however).

Contributing back bugfixes and improvements is polite and encouraged but not
required.

## Contact
Written by:
Jim McCann (jmccann@cs.cmu.edu) [http://tchow.com]
Additional text by Roger Dannenberg (rbd@cs.cmu.edu)
Feel free to email.
