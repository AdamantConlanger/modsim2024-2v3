from auxiliaries.simulate import simulate
from auxiliaries.cartesian_product import cartesian_product

model = input("which model?\n")

if model == "sfu":
    from sfu.single_run import perform_program
elif model == "stu":
    from stu.single_run import perform_program
elif model == "str":
    from str.single_run import perform_program
elif model == "stq":
    from stq.single_run import perform_program
elif model == "stp":
    print("model not available yet")
    exit()
    from stp.single_run import perform_program
elif model == "cfu":
    from cfu.single_run import perform_program
elif model == "ctu":
    from ctu.single_run import perform_program
elif model == "ctr":
    from ctr.single_run import perform_program
elif model == "ctq":
    from ctq.single_run import perform_program
elif model == "ctp":
    print("model not available yet")
    exit()
    from ctp.single_run import perform_program
else:
    print("no model found")
    exit()
perform_program(simulate, cartesian_product)
