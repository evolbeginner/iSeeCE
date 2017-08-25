class Dir
    def self.mkdirs(path)
        if(!File.directory?(path))
            if(!mkdirs(File.dirname(path)))
                return false;
            end
            mkdir(path)
        end
        return true
    end
end


################################################################################
def mkdir_with_force(outdir, is_force=false, is_tolerate=false)
  if outdir.class != String
    raise "outdir wrong? Exiting ......"
  end

  if ! Dir.exists?(outdir)
    `mkdir -p #{outdir}`
  else
    if is_tolerate
      ;
    elsif is_force
      `rm -rf #{outdir}`
      `mkdir -p #{outdir}`
    else
      raise "The outdir #{outdir} has already existed!"
    end
  end
end


