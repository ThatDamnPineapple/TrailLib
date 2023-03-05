using Microsoft.Xna.Framework;
using Microsoft.Xna.Framework.Graphics;
using System;
using System.Collections.Generic;
using Terraria;
using Terraria.Utilities;
using System.Linq;
using static TrailLib.PolygonUpscaling;

namespace TrailLib
{
    public class PrimitiveQuadrilateral : IDisposable
    {
        private readonly Primitives primitives;

        private const int verticesCount = 4;
        private const int indicesCount = 6;

        private readonly BasicEffect baseEffect;

        public Color color;

        /// <summary>
        /// Array of positions that define the trail. NOTE: Positions[Positions.Length - 1] is assumed to be the start (e.g. Projectile.Center) and Positions[0] is assumed to be the end.
        /// </summary>
        public Vector2[] Vertices {
            get => vertices;
            set {
                if (value.Length != verticesCount)
                {
                    throw new ArgumentException("Array of positions was a different length than the expected result!");
                }

                vertices = value;
            }
        }
        private Vector2[] vertices;

        public PrimitiveQuadrilateral(Color? color = null)
        {
            /* A---B
             * |  /|
             * C / D
             * Bozo alert
             */
            this.color = color ?? Color.White;
            primitives = new Primitives(Main.graphics.GraphicsDevice, verticesCount, indicesCount);
            baseEffect = new BasicEffect(Main.graphics.GraphicsDevice) { VertexColorEnabled = true, TextureEnabled = false };
        }

        private void GenerateMesh(out VertexPositionColorTexture[] vertices, out short[] indices)
        {
            Color col = Color.White;

            //HARDCODING IS THE BEST! SAY IT WITH ME! HARDCODING! RULES!
            vertices = new VertexPositionColorTexture[verticesCount] 
            {   new VertexPositionColorTexture(Vertices[0].Vec3(), color, Vector2.Zero), //Top left
                new VertexPositionColorTexture(Vertices[1].Vec3(), color, Vector2.UnitX), //Top right
                new VertexPositionColorTexture(Vertices[2].Vec3(), color, Vector2.UnitY), //Bottom left
                new VertexPositionColorTexture(Vertices[3].Vec3(), color, Vector2.One) //Bottom right
            };

            /* 0---1
             * |  /|
             * |/  |
             * 2---3
             */

            indices = new short[indicesCount] 
            {
                (short)1, (short)0,(short)2,
                (short)2, (short)3,(short)1,
            };
        }

        private void SetupMesh()
        {
            GenerateMesh(out VertexPositionColorTexture[] mainVertices, out short[] mainIndices);
            primitives.SetVertices(mainVertices);
            primitives.SetIndices(mainIndices);
        }

        public void Render(Effect effect = null, Vector2? offset = null)
        {
            Vector2 offset_ = offset.GetValueOrDefault();
            Render(effect, Matrix.CreateTranslation(offset_.Vec3()));
        }

        public void Render(Effect effect = null, Matrix? translation = null)
        {
            if (Vertices == null && !(primitives?.IsDisposed ?? true))
            {
                return;
            }

            Main.instance.GraphicsDevice.RasterizerState = RasterizerState.CullNone;
            SetupMesh();
            if (!translation.HasValue)
                translation = Matrix.CreateTranslation(0, 0, 0);

            if (effect == null)
                effect = baseEffect;

            primitives.Render(effect, translation.Value);
        }

        public void RenderWithView(Matrix view, Effect effect = null, Matrix? translation = null)
        {
            if (Vertices == null && !(primitives?.IsDisposed ?? true))
                return;

            Main.instance.GraphicsDevice.RasterizerState = RasterizerState.CullNone;
            SetupMesh();
            if (!translation.HasValue)
                translation = Matrix.CreateTranslation(0, 0, 0);

            if (effect == null)
                effect = baseEffect;

            primitives.Render(effect, translation.Value, view);
        }

        public void Dispose()
        {
            primitives?.Dispose();
        }
    }

    public class PrimitivePolygonOutline : IDisposable
    {
        private readonly Primitives primitives;

        private readonly int maxVertexCount;
        private readonly BasicEffect baseEffect;

        public Color color;
        public float outlineThickness;
        private Vector2 shapeCenter;

        public PolygonUpscalingAlgorithm outlineUpscaling;

        /// <summary>
        /// Array of positions that define the trail. NOTE: Positions[Positions.Length - 1] is assumed to be the start (e.g. Projectile.Center) and Positions[0] is assumed to be the end.
        /// </summary>
        public Vector2[] Vertices {
            get => vertices;
            set {
                if (value.Length != maxVertexCount)
                {
                    throw new ArgumentException("Array of positions was a different length than the expected result!");
                }

                vertices = value;


                Vector2 highestCoordinates = new Vector2(float.MinValue);
                Vector2 lowestCoordinates = new Vector2(float.MaxValue);
                foreach (Vector2 vertex in vertices)
                {
                    highestCoordinates.X = Math.Max(vertex.X, highestCoordinates.X);
                    highestCoordinates.Y = Math.Max(vertex.Y, highestCoordinates.Y);
                    lowestCoordinates.X = Math.Min(vertex.X, lowestCoordinates.X);
                    lowestCoordinates.Y = Math.Min(vertex.Y, lowestCoordinates.Y);
                }

                shapeCenter = (highestCoordinates + lowestCoordinates) / 2f;
            }
        }
        private Vector2[] vertices;

        /// <summary>
        /// THIS PROBABLY WONT WORK WITH CONVEX SHAPES
        /// DONT EVEN TRY I WILL NOT FIX THE UPSCALING
        /// </summary>
        /// <param name="pointCount"></param>
        /// <param name="color"></param>
        public PrimitivePolygonOutline(int pointCount, float outlineThickness, Color? color = null, PolygonUpscalingAlgorithm upscaler = null)
        {
            maxVertexCount = pointCount;
            this.color = color ?? Color.White;
            this.outlineThickness = outlineThickness;
            outlineUpscaling = upscaler ?? UpscalePolygonFromCenter;

            //https://media.discordapp.net/attachments/802291445360623686/1036075596888952883/unknown.png

            primitives = new Primitives(Main.graphics.GraphicsDevice, maxVertexCount * 2, maxVertexCount * 6);
            baseEffect = new BasicEffect(Main.graphics.GraphicsDevice) { VertexColorEnabled = true, TextureEnabled = false };
        }


        private void GenerateMesh(out VertexPositionColorTexture[] vertices, out short[] indices)
        {
            vertices = new VertexPositionColorTexture[maxVertexCount * 2];

            Vector2[] outwardsVertices = outlineUpscaling(Vertices.ToList(), outlineThickness).ToArray();

            // k = 0 indicates starting at the end of the trail (furthest from the origin of it).
            for (int k = 0; k < Vertices.Length; k++)
            {
                // 1 at k = Positions.Length - 1 (start) and 0 at k = 0 (end).
                float factorAlongOutline = (float)k / (Vertices.Length - 1);

                // Uses the trail width function to decide the width of the trail at this point (if no function, use 

                Vector2 inner = Vertices[k];
                Vector2 outer = outwardsVertices[k];

                //No nerd shit. Tldr X coordinate = how far along the outline from the first vertex to the last.
                //Y coordinate = is it on the outside or the inside of the outline
                Vector2 innerTexCoords = new Vector2(factorAlongOutline, 0);
                Vector2 outerTexCoords = new Vector2(factorAlongOutline, 1);

                //https://media.discordapp.net/attachments/802291445360623686/1036075596888952883/unknown.png

                vertices[k] = new VertexPositionColorTexture(outer.Vec3(), color, outerTexCoords);
                vertices[k + maxVertexCount] = new VertexPositionColorTexture(inner.Vec3(), color, innerTexCoords);
            }

            indices = new short[maxVertexCount * 6];

            /* Now, we have to loop through the indices to generate triangles.
             */
            for (short k = 0; k < maxVertexCount; k++)
            {
                short loopPoint = (short)(maxVertexCount * 2);
                //https://media.discordapp.net/attachments/802291445360623686/1036076440514482176/unknown.png

                short nextOuterIndex = (short)(k + 1);
                if (nextOuterIndex >= maxVertexCount)
                    nextOuterIndex = 0;



                indices[k * 6] = (short)(k % loopPoint);
                indices[k * 6 + 1] = (short)((k + maxVertexCount)); //First triangle
                indices[k * 6 + 2] = (short)(nextOuterIndex + maxVertexCount);

                indices[k * 6 + 3] = (short)(k);
                indices[k * 6 + 4] = nextOuterIndex; //Second triangle
                indices[k * 6 + 5] = (short)(nextOuterIndex + maxVertexCount);
            }
        }

        private void SetupMesh()
        {
            GenerateMesh(out VertexPositionColorTexture[] mainVertices, out short[] mainIndices);
            primitives.SetVertices(mainVertices);
            primitives.SetIndices(mainIndices);
        }

        public void Render(Effect effect = null, Vector2? offset = null)
        {
            Vector2 offset_ = offset.GetValueOrDefault();
            Render(effect, Matrix.CreateTranslation(offset_.Vec3()));
        }

        public void Render(Effect effect = null, Matrix? translation = null)
        {
            if (Vertices == null && !(primitives?.IsDisposed ?? true))
            {
                return;
            }

            Main.instance.GraphicsDevice.RasterizerState = RasterizerState.CullNone;
            SetupMesh();
            if (!translation.HasValue)
                translation = Matrix.CreateTranslation(0, 0, 0);

            if (effect == null)
                effect = baseEffect;

            primitives.Render(effect, translation.Value);
        }

        public void RenderWithView(Matrix view, Effect effect = null, Matrix? translation = null)
        {
            if (Vertices == null && !(primitives?.IsDisposed ?? true))
                return;

            Main.instance.GraphicsDevice.RasterizerState = RasterizerState.CullNone;
            SetupMesh();
            if (!translation.HasValue)
                translation = Matrix.CreateTranslation(0, 0, 0);

            if (effect == null)
                effect = baseEffect;

            primitives.Render(effect, translation.Value, view);
        }

        public void Dispose()
        {
            primitives?.Dispose();
        }
    }

    public class PrimitiveClosedLoop : IDisposable
    {
        private readonly Primitives primitives;

        private readonly int maxPointCount;

        private readonly TrailWidthFunction trailWidthFunction;

        private readonly TrailColorFunction trailColorFunction;

        private readonly BasicEffect baseEffect;

        /// <summary>
        /// Array of positions that define the trail. NOTE: Positions[Positions.Length - 1] is assumed to be the start (e.g. Projectile.Center) and Positions[0] is assumed to be the end.
        /// </summary>
        public Vector2[] Positions {
            get => positions;
            set {
                if (value.Length != maxPointCount)
                {
                    throw new ArgumentException("Array of positions was a different length than the expected result!");
                }

                positions = value;
            }
        }

        private Vector2[] positions;

        private const float defaultWidth = 16;

        public PrimitiveClosedLoop(int maxPointCount, TrailWidthFunction trailWidthFunction, TrailColorFunction trailColorFunction)
        {
            this.maxPointCount = maxPointCount;
            this.trailWidthFunction = trailWidthFunction;
            this.trailColorFunction = trailColorFunction;

            /* A---B---C
             * |  /|  /|
             * D / E / F
             * |/  |/  |
             * G---H---I
             * 
             * Let D, E, F, etc. be the set of n points that define the trail.
             * Since each point generates 2 vertices, there are 2n vertices, plus the tip's count.
             * 
             * As for indices - in the region between 2 defining points there are 2 triangles.
             * The amount of regions in the whole trail are given by n, so there are 2n triangles for n points.
             * Finally, since each triangle is defined by 3 indices, there are 6n  indices, plus the tip's count.
             */

            primitives = new Primitives(Main.graphics.GraphicsDevice, maxPointCount * 2 + 2, 6 * maxPointCount);
            baseEffect = new BasicEffect(Main.graphics.GraphicsDevice)
            {
                VertexColorEnabled = true,
                TextureEnabled = false
            };
        }

        private void GenerateMesh(out VertexPositionColorTexture[] vertices, out short[] indices)
        {
            VertexPositionColorTexture[] verticesTemp = new VertexPositionColorTexture[maxPointCount * 2 + 2];

            short[] indicesTemp = new short[maxPointCount * 6];

            //Generating vertices is the same as for prim trails
            // k = 0 indicates starting at the end of the trail (furthest from the origin of it).
            for (int k = 0; k <= Positions.Length; k++)
            {
                // 1 at k = Positions.Length - 1 (start) and 0 at k = 0 (end).
                float factorAlongTrail = (float)k / (Positions.Length);

                // Uses the trail width function to decide the width of the trail at this point (if no function, use 
                float width = trailWidthFunction?.Invoke(factorAlongTrail) ?? defaultWidth;

                Vector2 current = Positions[k % positions.Length];
                Vector2 next = Positions[(k + 1) % positions.Length];

                Vector2 normalToNext = (next - current).SafeNormalize(Vector2.Zero);
                Vector2 normalPerp = normalToNext.RotatedBy(MathHelper.PiOver2);

                /* A
                 * |
                 * B---D
                 * |
                 * C
                 * 
                 * Let B be the current point and D be the next one.
                 * A and C are calculated based on the perpendicular vector to the normal from B to D, scaled by the desired width calculated earlier.
                 */

                Vector2 a = current + (normalPerp * width);
                Vector2 c = current - (normalPerp * width);

                /* Texture coordinates are calculated such that the top-left is (0, 0) and the bottom-right is (1, 1).
                 * To achieve this, we consider the Y-coordinate of A to be 0 and that of C to be 1, while the X-coordinate is just the factor along the trail.
                 * This results in the point last in the trail having an X-coordinate of 0, and the first one having a Y-coordinate of 1.
                 */
                Vector2 texCoordA = new Vector2(factorAlongTrail, 0);
                Vector2 texCoordC = new Vector2(factorAlongTrail, 1);

                // Calculates the color for each vertex based on its texture coordinates. This acts like a very simple shader (for more complex effects you can use the actual shader).
                Color colorA = trailColorFunction?.Invoke(factorAlongTrail) ?? Color.White;
                Color colorC = trailColorFunction?.Invoke(factorAlongTrail) ?? Color.White;

                /* 0---1---2
                 * |  /|  /|
                 * A / B / C
                 * |/  |/  |
                 * 3---4---5
                 * 
                 * Assuming we want vertices to be indexed in this format, where A, B, C, etc. are defining points and numbers are indices of mesh points:
                 * For a given point that is k positions along the chain, we want to find its indices.
                 * These indices are given by k for the above point and k + n + 1 for the below point. (+1 cuz we have the doubled final points
                 */

                verticesTemp[k] = new VertexPositionColorTexture(a.Vec3(), colorA, texCoordA);
                verticesTemp[k + maxPointCount + 1] = new VertexPositionColorTexture(c.Vec3(), colorC, texCoordC);
            }

            /* Now, we have to loop through the indices to generate triangles.
             * Looping to maxPointCount brings us halfway to the end; it covers the top row (excluding the last point on the top row).
             */
            for (short k = 0; k < maxPointCount - 1; k++)
            {
                /* 0---1
                 * |  /|
                 * A / B
                 * |/  |
                 * 2---3
                 * 
                 * This illustration is the most basic set of points (where n = 2).
                 * In this, we want to make triangles (2, 3, 1) and (1, 0, 2).
                 * Generalising this, if we consider A to be k = 0 and B to be k = 1, then the indices we want are going to be (k + n, k + n + 1, k + 1) and (k + 1, k, k + n)
                 */

                indicesTemp[k * 6] = (short)(k + maxPointCount + 1);
                indicesTemp[k * 6 + 1] = (short)(k + maxPointCount + 2);
                indicesTemp[k * 6 + 2] = (short)(k + 1);
                indicesTemp[k * 6 + 3] = (short)(k + 1);
                indicesTemp[k * 6 + 4] = k;
                indicesTemp[k * 6 + 5] = (short)(k + maxPointCount + 1);
            }

            //Set the final triangles to loop the strip
            indicesTemp[(maxPointCount - 1) * 6] = (short)(maxPointCount * 2); //maxPointCount - 1 + maxPointCount + 1
            indicesTemp[(maxPointCount - 1) * 6 + 1] = (short)(maxPointCount * 2 + 1); //maxPointCount - 1 + maxPointCount + 2
            indicesTemp[(maxPointCount - 1) * 6 + 2] = (short)0;
            indicesTemp[(maxPointCount - 1) * 6 + 3] = (short)0;
            indicesTemp[(maxPointCount - 1) * 6 + 4] = (short)(maxPointCount - 1);
            indicesTemp[(maxPointCount - 1) * 6 + 5] = (short)(maxPointCount * 2);//maxPointCount - 1 + maxPointCount + 1

            // The next available index will be the next value after the count of points (starting at 0).
            vertices = verticesTemp;
            // Maybe we could use an array instead of a list for the indices, if someone figures out how to add indices to an array properly.
            indices = indicesTemp;
        }

        private void SetupMeshes()
        {
            GenerateMesh(out VertexPositionColorTexture[] mainVertices, out short[] mainIndices);
            primitives.SetVertices(mainVertices);
            primitives.SetIndices(mainIndices);
        }

        public void Render(Effect effect = null, Vector2? offset = null)
        {
            Vector2 offset_ = offset.GetValueOrDefault();
            Render(effect, Matrix.CreateTranslation(offset_.Vec3()));
        }

        public void Render(Effect effect = null, Matrix? translation = null)
        {
            if (Positions == null && !(primitives?.IsDisposed ?? true))
            {
                return;
            }

            Main.instance.GraphicsDevice.RasterizerState = RasterizerState.CullNone;

            SetupMeshes();
            if (!translation.HasValue)
                translation = Matrix.CreateTranslation(-Main.screenPosition.Vec3());

            if (effect == null)
            {
                effect = baseEffect;
            }

            primitives.Render(effect, translation.Value);
        }

        public void Dispose()
        {
            primitives?.Dispose();
        }

        public void SetPositionsCircle(Vector2 Center, float radius, float rotationOffset = 0f)
        {
            Vector2[] positionsTemp = new Vector2[maxPointCount];

            for (int i = 0; i < maxPointCount; i++)
            {
                float rotation = i / (float)maxPointCount * MathHelper.TwoPi + rotationOffset;
                positionsTemp[i] = Center + rotation.ToRotationVector2() * radius;
            }

            Positions = positionsTemp;
        }
    }

    public static class PolygonUpscaling
    {
        public delegate List<Vector2> PolygonUpscalingAlgorithm(List<Vector2> vertices, float sizeUp);

        public static List<Vector2> UpscalePolygonFromCenter(List<Vector2> vertices, float sizeUp)
        {
            List<Vector2> upscaledVertices = new List<Vector2>();

            Vector2 highestCoordinates = new Vector2(float.MinValue);
            Vector2 lowestCoordinates = new Vector2(float.MaxValue);

            foreach (Vector2 vertex in vertices)
            {
                highestCoordinates.X = Math.Max(vertex.X, highestCoordinates.X);
                highestCoordinates.Y = Math.Max(vertex.Y, highestCoordinates.Y);
                lowestCoordinates.X = Math.Min(vertex.X, lowestCoordinates.X);
                lowestCoordinates.Y = Math.Min(vertex.Y, lowestCoordinates.Y);
            }

            Vector2 center = (highestCoordinates + lowestCoordinates) / 2f;

            for (int i = 0; i < vertices.Count; i++)
            {
                Vector2 fromCenter = (vertices[i] - center);
                upscaledVertices.Add(center + fromCenter.SafeNormalize(Vector2.One) * (fromCenter.Length() + sizeUp));
            }
            return upscaledVertices;
        }

        public static List<Vector2> UpscalePolygonByCombiningPerpendiculars(List<Vector2> vertices, float sizeUp)
        {
            List<Vector2> upscaledVertices = new List<Vector2>();

            Vector2 previousPerpendicular = (vertices[vertices.Count - 1] - vertices[0]).RotatedBy(MathHelper.PiOver2).SafeNormalize(Vector2.Zero);

            for (int i = 0; i < vertices.Count; i++)
            {
                Vector2 nextPerpendicular = (vertices[i] - vertices[(i + 1) % vertices.Count]).RotatedBy(MathHelper.PiOver2).SafeNormalize(Vector2.Zero);
                upscaledVertices.Add(vertices[i] + (nextPerpendicular + previousPerpendicular) * sizeUp);
                previousPerpendicular = nextPerpendicular;
            }

            return upscaledVertices;
        }

        public static T Next<T>(this UnifiedRandom random, IEnumerable<T> collection)
        {
            return collection.ElementAt(random.Next(collection.Count()));
        }

        public static List<Vector2> UpscalePolygonByCombininParallels(List<Vector2> vertices, float sizeUp)
        {
            List<Vector2> upscaledVertices = new List<Vector2>();


            Vector2 previousParallel = (vertices[vertices.Count - 1] - vertices[0]).SafeNormalize(Vector2.Zero);

            for (int i = 0; i < vertices.Count; i++)
            {
                Vector2 nextParallel = (vertices[i] - vertices[(i + 1) % vertices.Count]).SafeNormalize(Vector2.Zero);
                upscaledVertices.Add(vertices[i] + (nextParallel + previousParallel) * sizeUp);
                previousParallel = nextParallel;
            }

            return upscaledVertices;
        }
    }
}
