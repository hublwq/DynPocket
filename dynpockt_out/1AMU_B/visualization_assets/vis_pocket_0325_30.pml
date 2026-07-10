load ../conformations_representative/rep_frame_0325.pdb, frame_0325
hide all
show cartoon, frame_0325
color white, frame_0325
select target_pocket, frame_0325 and resi 255+170+171+172+301+299+303+300+210+211+212+307+213+279+278+308+276+216
show sticks, target_pocket
util.cbay target_pocket
show surface, target_pocket
set transparency, 0.5, target_pocket
zoom target_pocket
