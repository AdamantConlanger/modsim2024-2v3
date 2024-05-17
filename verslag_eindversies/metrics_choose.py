from auxiliaries.simulate import simulate
from auxiliaries.cartesian_product import cartesian_product

model = input("which model?\n")

if model == "sfu":
    print("metrics not possible for model")
    exit()
    from sfu.metrics import perform_program
elif model == "stu":
    print("model not available yet")
    exit()
    from stu.metrics import perform_program
elif model == "str":
    print("model not available yet")
    exit()
    from str.metrics import perform_program
elif model == "stq":
    from stq.metrics import perform_program
elif model == "stp":
    print("model not available yet")
    exit()
    from stp.metrics import perform_program
elif model == "efu":
    print("metrics not possible for model")
    exit()
    from cfu.metrics import perform_program
elif model == "ctu":
    print("model not available yet")
    exit()
    from ctu.metrics import perform_program
elif model == "ctr":
    print("model not available yet")
    exit()
    from ctr.metrics import perform_program
elif model == "ctq":
    from ctq.metrics import perform_program
elif model == "ctp":
    print("model not available yet")
    exit()
    from ctp.metrics import perform_program
else:
    print("no model found")
    exit()
perform_program(simulate, cartesian_product)
