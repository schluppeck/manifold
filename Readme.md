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


matlab code
-----------

a couple of simple matlab functions for loading in vtk files and basics of displaying.

* use `loadSurfVTK()` to get meshes into matlab
* `patch()` to render, e.g.

	fname = 'rh.pial.vtk'; % white matter surface
	s = loadSurfVTK(fname);
    figure, patch('vertices',s.vtcs, 'faces',s.tris, ...
     'facevertexcdata', ones(s.Nvtcs,1)), colormap(jet)
    shading interp
	material dull
	light
	axis equal, axis vis3d, axis off

