from auxiliaries.simulate import simulate
from auxiliaries.cartesian_product import cartesian_product

model = input("which model?\n")

if model == "sfu":
    from sfu.multiple_runs import perform_program
elif model == "stu":
    from stu.multiple_runs import perform_program
elif model == "str":
    from str.multiple_runs import perform_program
elif model == "stq":
    from stq.multiple_runs import perform_program
elif model == "stp":
    print("model not available yet")
    exit()
    from stp.multiple_runs import perform_program
elif model == "cfu":
    from cfu.multiple_runs import perform_program
elif model == "ctu":
    from ctu.multiple_runs import perform_program
elif model == "ctr":
    from ctr.multiple_runs import perform_program
elif model == "ctq":
    from ctq.multiple_runs import perform_program
elif model == "ctp":
    print("model not available yet")
    exit()
    from ctp.multiple_runs import perform_program
else:
    print("no model found")
    exit()
perform_program(simulate, cartesian_product)
