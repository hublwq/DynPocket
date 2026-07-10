load ../conformations_representative/rep_frame_0000.pdb, ref_struct
bg_color white
hide all
show cartoon, ref_struct
color gray80, ref_struct
select core_pocket, ref_struct and resi 308+303+359+216+304+386+212+255+301+299+101+307+207+306+177+279+305+300+278+358
show surface, core_pocket
color red, core_pocket
set transparency, 0.4, core_pocket
zoom core_pocket
