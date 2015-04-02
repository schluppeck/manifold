manifold
========

general
-------

matlab code, info and some snippets of code for looking at functional and anatomical MRI data on the cortical surface with colleagues at the UoN School of Mathematics.


data
----

anatomical data derived from freesurfer segementations of a sample subject.

to get vtk file for e.g. the right hemisphere, you can `cd` to the subject's directory `${SUBJECTS_DIR}/xx/surf/` and convert all the freesurfer format files for **lh** and **rh** by looping	

	# make sure extglob option in BASH is on
	# extended globbing of files...
	shopt -s extglob 

    for fname in +(rh|lh).+([a-z]) ; do echo "convert $fname to $fname.vtk"; done

    for fname in +(rh|lh).+([a-z]) ; do mris_convert $fname $fname.vtk; done

note that for scalar data for which there is only one value per vertex is (eg. curv, thickness or area files), `mris_convert` needs some additional information. in this case you need to specify the surface from which data was derived with an additional option.

	# add -c option
	mris_convert -c lh.curv lh.white lh.curv.vtk
	mris_convert -c rh.curv rh.white lh.rurv.vtk
	

matlab code
-----------

a couple of simple matlab functions for loading in vtk files and basics of displaying.

inspect `dataInfoS1` to get a sense of how to use various functions. but here as a quick example:

    fname = 'rh.pial.vtk'; % white matter surface
    s = loadSurfVTK(fname);
    figure
    renderSurf(s) % wrapper around |patch|
    colormap(jet)
    shading interp
    material dull
    light
    axis equal, axis vis3d, axis off

