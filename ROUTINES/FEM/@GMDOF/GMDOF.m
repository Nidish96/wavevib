    classdef GMDOF
  %GMDOF creates a generic mdof model class.
  properties
    Nds		% Nodal locations (only contacting dofs) [x y z...]
    Nn		% Number of nodes

    dpn		% Number of dofs per node

    Qm 		% node-nlpoint interpolation matrix
    Tm		% nlpoints-nodes integration matrix

    fcont 	% Contact Function
    z           % Slider States for fcont
  end

  methods
    function m = GMDOF(nds, dpn, Qm, Tm)
    %GMDOF initializes the GMDOF object.
    %
    % USAGE:
    %  m = GMDOF(Nds, Tri, Quad);
      % Nodes & Dofs per node
      m.Nds = nds;
      m.Nn  = size(nds, 1);
      m.dpn = dpn;

      m.Qm = Qm;
      m.Tm = Tm;
    end

    function m = SETCFUN(m, fcont, varargin)
      m.fcont = fcont;
      if length(varargin)==1
          m.z = varargin{1};
      end
    end

    function m = SINGLEPREC(m)
      m.Nds    = single(m.Nds);
      m.Nn     = uint16(m.Nn);
      m.dpn = uint16(m.dpn);
      
      m.Qm = m.Qm;
      m.Tm = m.Tm;
      
      m.z = single(m.z);
    end
  end
end
