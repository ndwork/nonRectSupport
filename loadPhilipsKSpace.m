
function kData = loadPhilipsKSpace( dataFile )
  % Taken from example_loadmixes.m

  kspace = loadRawKspace( dataFile );

  nCoils = kspace.kspace_properties.number_of_coil_channels;
  nDynamics = kspace.kspace_properties.number_of_dynamic_scans;
  nMixes = kspace.kspace_properties.number_of_mixes;

  kData = zeros( kspace.kspace_properties.kx_range(2) - kspace.kspace_properties.kx_range(1)+1,... kx size
    max(kspace.ky) - min(kspace.ky) + 1, ...
    max(kspace.kz) - min(kspace.kz) + 1, ...
    nCoils, nMixes, nDynamics );

  selection_proper_lines = strcmp( kspace.typ, 'STD' );

  for index = find(selection_proper_lines)'
    kData(  :,... fill entire kx (readout dir.)
      kspace.ky(index) - min(kspace.ky) + 1,...
      kspace.kz(index) - min(kspace.kz) + 1,...
      kspace.chan(index) + 1, ...
      kspace.mix(index) + 1, ...
      kspace.dyn(index) + 1) ...
      = kspace.complexdata{ index };
    if rem(index,round(sum(selection_proper_lines)/40)) == 0
      disp(['Loading Progress: ' num2str(round((index/sum(selection_proper_lines))*100)) ' %']);
    end
  end

end
