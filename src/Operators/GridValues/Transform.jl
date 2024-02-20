function (tf::Transform{T})(n::Integer) where T <: GridValues
    grid = tf.basis.GD.D.map.(tf.basis.GD.grid(n))
    plan_transform(tf.basis,n,grid,x -> x)
end