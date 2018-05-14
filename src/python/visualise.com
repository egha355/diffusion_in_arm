#@exnodes=<./Laplace.part*.exnode>;
#@exelems=<./Laplace.part*.exelem>;
#foreach $filename (@exnodes) {
#    print "Reading $filename\n";
#    gfx read node "$filename";
#}
#foreach $filename (@exelems) {
#    print "Reading $filename\n";
#    gfx read elem "$filename";
#}

for ($i=0;$i<810;$i=$i+100) 
#Read in the sequence of nodal positions.
  {
     $filename = sprintf("./upperlimb2material/MainTime_%01d.part0.exnode", $i);
     $time = $i
     print "Reading $filename time $time\n";
     gfx read node "$filename" time $i;
  }

#Read in the element description
gfx read element ./upperlimb2material/MainTime_0.part0.exelem;

#gfx read node Laplace.part0.exnode
#gfx read elem Laplace.part0.exelem

gfx define faces egroup DiffusionRegion
gfx modify window 1 view interest_point 1.0,0.5,0.0 eye_point 1.0,0.5,5.0 up_vector 0.0,1.0,0.0
gfx modify spectrum default clear overwrite_colour
gfx modify spectrum default linear reverse range 0.0 100.0 extend_above extend_below rainbow colour_range 0 1 component 1
#gfx modify spectrum default linear reverse range 0.0 1.0 extend_above extend_below banded number_of_bands 20 band_ratio 0.06 component 1
#gfx modify g_element DiffusionRegion surfaces select_on material default data U spectrum default selected_material default_selected render_shaded;
#gfx modify g_element DiffusionRegion node_points glyph sphere general size "0.05*0.05*0.05" centre 0,0,0  label U
gfx modify g_element DiffusionRegion general clear circle_discretization 6 default_coordinate Coordinate element_discretization "1*1*1" native_discretization none;
gfx modify g_element DiffusionRegion node_points glyph sphere general size "0.08*0.08*0.08" centre 0,0,0 font default select_on material default data U spectrum default selected_material default_selected;
#gfx modify spectrum default clear overwrite_colour;
gfx modify g_element DiffusionRegion lines select_on material default selected_material default_selected
gfx modify g_element DiffusionRegion surfaces select_on material default data U spectrum default selected_material default_selected render_shaded

gfx draw axes
gfx edit scene
gfx create window 1

gfx cre mat copper ambient 1 0.2 0 diffuse 0.6 0.3 0 specular 0.7 0.7 0.5 shininess 0.3;
gfx create colour_bar label_material black spectrum default material copper;
gfx modify g_element "/" point glyph colour_bar general size "1*1*1" centre 0,0,0 select_on material copper selected_material copper normalised_window_fit_left;
