BEGIN{INSERT_NTIMEPOINTSM_NTIMEVARSM=1}
{

  if($4=="(NTIMEPOINTSM)[1]:")
    {
      INSERT_NTIMEPOINTSM_NTIMEVARSM=0;
    }
  else if($4=="(NTIMEVARSM)[1]:")
    {
      INSERT_NTIMEPOINTSM_NTIMEVARSM=0;
    }
   else if($7=="(NYOM)[1]:" && INSERT_NTIMEPOINTSM_NTIMEVARSM)
    {
      print(" Max# time points          (NTIMEPOINTSM)[1]:         1");
      print(" Max# time variables         (NTIMEVARSM)[1]:         1");
    }

  print $0;
}
END{}
