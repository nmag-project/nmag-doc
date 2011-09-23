import os

def test_all():
  from nsim.testtools import get_directory_of, run_nsim_in_dir
  exit_status = run_nsim_in_dir(get_directory_of(__file__), "run.py")
  assert exit_status == 0

if __name__ == "__main__":
  test_all()

