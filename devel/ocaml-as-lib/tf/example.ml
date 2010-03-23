
let demo () = 
  (* We link in Unix as an example for some external module *)
  let s = Unix.stat "/etc/passwd" in
  Printf.printf "The size of /etc/passed is %d bytes.\n%!" s.Unix.st_size
;;

Callback.register "caml_callback_demo" demo;;
