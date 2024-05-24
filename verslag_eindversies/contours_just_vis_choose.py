from auxiliaries.simulate import simulate
from auxiliaries.cartesian_product import cartesian_product

model = input("which model?\n")

if model == "sfu":
    print("contours not possible for model")
    exit()
    from sfu.contours_just_vis import perform_program
elif model == "stu":
    print("model not available yet")
    exit()
    from stu.contours_just_vis import perform_program
elif model == "str":
    from str.contours_just_vis import perform_program
elif model == "stq":
    print("model not available yet")
    exit()
    from stq.contours_just_vis import perform_program
elif model == "stp":
    print("model not available yet")
    exit()
    from stp.contours_just_vis import perform_program
elif model == "efu":
    print("contours not possible for model")
    exit()
    from cfu.contours_just_vis import perform_program
elif model == "ctu":
    print("model not available yet")
    exit()
    from ctu.contours_just_vis import perform_program
elif model == "ctr":
    from ctr.contours_just_vis import perform_program
elif model == "ctq":
    print("model not available yet")
    exit()
    from ctq.contours_just_vis import perform_program
elif model == "ctp":
    print("model not available yet")
    exit()
    from ctp.contours_just_vis import perform_program
else:
    print("no model found")
    exit()
perform_program(simulate, cartesian_product)
